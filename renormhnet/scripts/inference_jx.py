import os
import time
import subprocess as sp
import json

import argparse

from renormhnet import get_inference_config, rotate_back, abacus_parse
from renormhnet.preprocess import openmx_parse_overlap, get_rc, openmx_parse_overlap_jx, get_rcj #added by DHK
from renormhnet.inference import predict, predict_with_grad
from renormhnet.inference import plotjx

# for jx mode: task = [2,3,4] no overlap and eigenvalue evaluation
# 1 in list means to parse the overlap, 
# 2 means to get the local coordinate, 
# 3 means to predict the Hamiltonian, 
# 4 means to rotate the Hamiltonian back,
# 5 means to perform the sparse calculation. (jx mode: plot J vs. distance -> jx_plot.png)

def main():
    parser = argparse.ArgumentParser(description='Deep Hamiltonian')
    parser.add_argument('--config', default=[], nargs='+', type=str, metavar='N')
    args = parser.parse_args()

    print(f'User config name: {args.config}')
    config = get_inference_config(args.config)

    work_dir = os.path.abspath(config.get('basic', 'work_dir'))
    OLP_dir = os.path.abspath(config.get('basic', 'OLP_dir')) #overlap path
    interface = config.get('basic', 'interface')  # openmx | abacus
    abacus_suffix = str(config.get('basic', 'abacus_suffix', fallback='ABACUS'))
    task = json.loads(config.get('basic', 'task'))
    assert isinstance(task, list)
    eigen_solver = config.get('basic', 'eigen_solver')
    disable_cuda = config.getboolean('basic', 'disable_cuda')
    device = config.get('basic', 'device')
    huge_structure = config.getboolean('basic', 'huge_structure')
    restore_blocks_py = config.getboolean('basic', 'restore_blocks_py')
    gen_rc_idx = config.getboolean('basic', 'gen_rc_idx')
    gen_rc_by_idx = config.get('basic', 'gen_rc_by_idx')
    with_grad = config.getboolean('basic', 'with_grad')
    julia_interpreter = config.get('interpreter', 'julia_interpreter', fallback='')
    python_interpreter = config.get('interpreter', 'python_interpreter', fallback='')
    radius = config.getfloat('graph', 'radius')

    if 5 in task and interface not in ['jx', 'jxmnte']: # jx mode uses task 5 for plotting, no eigen solver
        if eigen_solver in ['sparse_jl', 'dense_jl']:
            assert julia_interpreter, "Please specify julia_interpreter to use Julia code to calculate eigenpairs"
        elif eigen_solver in ['dense_py']:
            assert python_interpreter, "Please specify python_interpreter to use Python code to calculate eigenpairs"
        else:
            raise ValueError(f"Unknown eigen_solver: {eigen_solver}")
    if 3 in task and not restore_blocks_py:
        assert julia_interpreter, "Please specify julia_interpreter to use Julia code to rearrange matrix blocks"

    if with_grad: # default False
        assert restore_blocks_py is True
        assert 4 not in task
        assert 5 not in task

    os.makedirs(work_dir, exist_ok=True)
    config.write(open(os.path.join(work_dir, 'config.ini'), "w"))


    if not restore_blocks_py: #default True
        cmd3_post = f"{julia_interpreter} " \
                    f"{os.path.join(os.path.dirname(os.path.dirname(__file__)), 'inference', 'restore_blocks.jl')} " \
                    f"--input_dir {work_dir} --output_dir {work_dir}"

    #no need to care about eigen solver 
    if eigen_solver == 'sparse_jl':
        cmd5 = f"{julia_interpreter} " \
               f"{os.path.join(os.path.dirname(os.path.dirname(__file__)), 'inference', 'sparse_calc.jl')} " \
               f"--input_dir {work_dir} --output_dir {work_dir} --config {config.get('basic', 'sparse_calc_config')}"
    elif eigen_solver == 'dense_jl':
        cmd5 = f"{julia_interpreter} " \
               f"{os.path.join(os.path.dirname(os.path.dirname(__file__)), 'inference', 'dense_calc.jl')} " \
               f"--input_dir {work_dir} --output_dir {work_dir} --config {config.get('basic', 'sparse_calc_config')}"
    elif eigen_solver == 'dense_py':
        cmd5 = f"{python_interpreter} " \
               f"{os.path.join(os.path.dirname(os.path.dirname(__file__)), 'inference', 'dense_calc.py')} " \
               f"--input_dir {work_dir} --output_dir {work_dir} --config {config.get('basic', 'sparse_calc_config')}"
    else:
        raise ValueError(f"Unknown eigen_solver: {eigen_solver}")

    print(f"\n~~~~~~~ 1.parse_Overlap\n")
    print(f"\n~~~~~~~ 2.get_local_coordinate\n")
    print(f"\n~~~~~~~ 3.get_pred_Hamiltonian\n")
    if not restore_blocks_py:
        print(f"\n~~~~~~~ 3_post.restore_blocks, command: \n{cmd3_post}\n")
    print(f"\n~~~~~~~ 4.rotate_back\n")
    if interface in ['jx', 'jxmnte']:
        print(f"\n~~~~~~~ 5.plot_jx\n")
    else:
        print(f"\n~~~~~~~ 5.sparse_calc, command: \n{cmd5}\n")

    # following files are required:
    # hdf5 file - overlaps.h5 - to find neighbors, which means local environment.
    # lat.dat
    # rlat.dat
    # site_positions.dat
    # orbital_types.dat  - but this should be changed as Heisenberg spin basis i.e. 0 for all elements
    # element.dat

    if 1 in task:
        begin = time.time()
        print(f"\n####### Begin 1.parse_Overlap")
        if interface == 'openmx':
            assert os.path.exists(os.path.join(OLP_dir, 'openmx.out')), "Necessary files could not be found in OLP_dir"
            assert os.path.exists(os.path.join(OLP_dir, 'output')), "Necessary files could not be found in OLP_dir"
            openmx_parse_overlap(OLP_dir, work_dir)
        elif interface == 'jx':
            assert os.path.exists(os.path.join(OLP_dir, 'openmx.out')), "Necessary files could not be found in OLP_dir"
            assert os.path.exists(os.path.join(OLP_dir, 'output')), "Necessary files could not be found in OLP_dir"
            #this part can be displaced for any bonds.
            openmx_parse_overlap_jx(OLP_dir, work_dir)    
        #elif interface == 'jxmnte':
        #    assert os.path.exists(os.path.join(OLP_dir, 'openmx.out')), "Necessary files could not be found in OLP_dir"
        #    assert os.path.exists(os.path.join(OLP_dir, 'output')), "Necessary files could not be found in OLP_dir"
        #    #this part can be displaced for any bonds.
        #    openmx_parse_overlap_jxmnte(OLP_dir, work_dir)     

            
        elif interface == 'abacus':
            print("Output subdirectories:", "OUT." + abacus_suffix)
            assert os.path.exists(os.path.join(OLP_dir, 'SR.csr')), "Necessary files could not be found in OLP_dir"
            assert os.path.exists(os.path.join(OLP_dir, f'OUT.{abacus_suffix}')), "Necessary files could not be found in OLP_dir"
            abacus_parse(OLP_dir, work_dir, data_name=f'OUT.{abacus_suffix}', only_S=True)
        assert os.path.exists(os.path.join(work_dir, "overlaps.h5"))
        assert os.path.exists(os.path.join(work_dir, "lat.dat"))
        assert os.path.exists(os.path.join(work_dir, "rlat.dat"))
        assert os.path.exists(os.path.join(work_dir, "site_positions.dat"))
        assert os.path.exists(os.path.join(work_dir, "orbital_types.dat"))
        assert os.path.exists(os.path.join(work_dir, "element.dat"))
        print('\n******* Finish 1.parse_Overlap, cost %d seconds\n' % (time.time() - begin))



    if not with_grad and 2 in task:
        begin = time.time()
        print(f"\n####### Begin 2.get_local_coordinate")
        #def get_rc(input_dir, output_dir, radius, r2_rand=False, gen_rc_idx=False, gen_rc_by_idx="", create_from_DFT=True, neighbour_file='overlaps.h5', if_require_grad=False, cart_coords=None):
        #if interface == 'jxmnte':
        if interface in ['jxmnte','jx']:
            def get_atoms(elementfile=os.path.join(work_dir, "element.dat")):
                with open(elementfile, "r") as f:
                    n = sum(1 for _ in f)
                line = ""
                for i in range(1,n+1):
                    line += f"{str(i)} "
                return line
            atomlabels = get_atoms()
            print("radius:", radius)
            cmdjx = f"{julia_interpreter} " \
                      f"{os.path.join(os.path.dirname(os.path.dirname(__file__)), 'preprocess', 'mk_hammat.jl')} " \
                      f"--atom {atomlabels} " \
                      f"--radius {radius} " \
                      f"--outputname  {work_dir}/Spin_Matrix.csv"
            print(cmdjx)
            sp.run(cmdjx, shell=True, capture_output=True, encoding="utf-8")                      
            cmd2jx = f"{julia_interpreter} " \
                      f"{os.path.join(os.path.dirname(os.path.dirname(__file__)), 'preprocess', 'jx_to_hdf5.jl')} " \
                      f"--input_dir {work_dir} --output_dir {work_dir}" 
            sp.run(cmd2jx, shell=True, capture_output=True, encoding="utf-8")                      
            
            get_rcj(work_dir, work_dir, radius=radius, gen_rc_idx=gen_rc_idx, gen_rc_by_idx=gen_rc_by_idx,
               create_from_DFT=config.getboolean('graph', 'create_from_DFT'), neighbour_file='Jhamiltonians.h5')
            os.remove(os.path.join(work_dir, "Jhamiltonians.h5"))
            os.remove(os.path.join(work_dir, "Spin_Matrix.csv"))
        else:
            get_rc(work_dir, work_dir, radius=radius, gen_rc_idx=gen_rc_idx, gen_rc_by_idx=gen_rc_by_idx,
               create_from_DFT=config.getboolean('graph', 'create_from_DFT'))               
        assert os.path.exists(os.path.join(work_dir, "rc.h5"))
        print('\n******* Finish 2.get_local_coordinate, cost %d seconds\n' % (time.time() - begin))

    if 3 in task:
        begin = time.time()
        print(f"\n####### Begin 3.get_pred_Hamiltonian")
        trained_model_dir = config.get('basic', 'trained_model_dir')
        if trained_model_dir[0] == '[' and trained_model_dir[-1] == ']':
            trained_model_dir = json.loads(trained_model_dir)
        if with_grad:
            predict_with_grad(input_dir=work_dir, output_dir=work_dir, disable_cuda=disable_cuda, device=device,
                              huge_structure=huge_structure, trained_model_dirs=trained_model_dir)
        else:
            predict(input_dir=work_dir, output_dir=work_dir, disable_cuda=disable_cuda, device=device,
                    huge_structure=huge_structure, restore_blocks_py=restore_blocks_py,
                    trained_model_dirs=trained_model_dir)
        if restore_blocks_py:
            if with_grad:
                assert os.path.exists(os.path.join(work_dir, "hamiltonians_grad_pred.h5"))
                assert os.path.exists(os.path.join(work_dir, "hamiltonians_pred.h5"))
            else:
                assert os.path.exists(os.path.join(work_dir, "rh_pred.h5"))
        else:
            capture_output = sp.run(cmd3_post, shell=True, capture_output=False, encoding="utf-8")
            assert capture_output.returncode == 0
            assert os.path.exists(os.path.join(work_dir, "rh_pred.h5"))
        print('\n******* Finish 3.get_pred_Hamiltonian, cost %d seconds\n' % (time.time() - begin))

    if 4 in task:
        begin = time.time()
        print(f"\n####### Begin 4.rotate_back")
        rotate_back(input_dir=work_dir, output_dir=work_dir)
        assert os.path.exists(os.path.join(work_dir, "hamiltonians_pred.h5"))
        print('\n******* Finish 4.rotate_back, cost %d seconds\n' % (time.time() - begin))
        
        
        cmd_hdf5_to_csv = f"{julia_interpreter} " \
            f"{os.path.join(os.path.dirname(os.path.dirname(__file__)), 'inference', 'hdf5_to_csv.jl')} " \
            f"--input {os.path.join(work_dir, 'hamiltonians_pred.h5')} " \
            f"--output {os.path.join(work_dir, 'Hamiltonian.csv')} " 
            
        sp.run(cmd_hdf5_to_csv, shell=True, capture_output=True, encoding="utf-8")

    if 5 in task:
        begin = time.time()
        if interface in ['jx', 'jxmnte']:
            print(f"\n####### Begin 5.plot_jx")
            assert os.path.exists(os.path.join(work_dir, "Hamiltonian.csv")), \
                "Hamiltonian.csv not found in work_dir; run task 4 first (hdf5_to_csv.jl)"
            cwd = os.getcwd()
            os.chdir(work_dir)  # plotjx() reads lat.dat, site_positions.dat, Hamiltonian.csv from cwd
            try:
                plotjx()
            finally:
                os.chdir(cwd)
            assert os.path.exists(os.path.join(work_dir, "jx_plot.png"))
            print('\n******* Finish 5.plot_jx, cost %d seconds\n' % (time.time() - begin))
        else:
            print(f"\n####### Begin 5.sparse_calc")
            capture_output = sp.run(cmd5, shell=True, capture_output=False, encoding="utf-8")
            assert capture_output.returncode == 0
            print('\n******* Finish 5.sparse_calc, cost %d seconds\n' % (time.time() - begin))

if __name__ == '__main__':
    main()
