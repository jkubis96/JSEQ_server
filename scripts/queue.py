turn = 'ON'
import os
import time
import sys
import subprocess

task_list_path=sys.argv[1]


def queueu(task_list_path, turn:str() = 'ON'):
    while(turn == 'ON'):
        try:
            with open(task_list_path) as f:
                tasks = f.read().splitlines()
            
            for task in tasks:
                
                #bashCommand = 'nohup ./docker_init projects/' + task + '/config &'

                command1 = subprocess.Popen(['nohup', './docker_init', 'projects/' + task + '/config', '&'])
                command1.wait()
                
                #os.system(bashCommand)
                try:
                    with open(task_list_path, 'r') as fin:
                        data = fin.read().splitlines()
                    with open(task_list_path, 'w') as final:
                        final.writelines('\n'.join(data[1:]))
                except:
                    print("Lack of task file")
                    time.sleep(60)

                    
        except:
            print("Lack of task file")
            time.sleep(60)
                
        
    
queueu(task_list_path, turn)


