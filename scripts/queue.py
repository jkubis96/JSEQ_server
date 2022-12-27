turn = 'ON'
import os
import time
import sys
task_list_path=sys.argv[1]


def queueu(task_list_path, turn:str() = 'ON'):
    while(turn == 'ON'):
        try:
            with open(task_list_path) as f:
                tasks = f.read().splitlines()
            
            for task in tasks:
                bashCommand = 'nohup ./docker_init $' + task + '&'
                os.system(bashCommand)
                tasks = tasks[tasks not in task]
                try:
                    with open(task_list_path, 'r') as fin:
                        data = fin.read().splitlines(True)
                    with open(task_list_path, 'w') as f:
                        f.writelines(data[1:])
                except:
                    print("Lack of task file")
                    time.sleep(60)

                    
        except:
            print("Lack of task file")
            time.sleep(60)
                
        
    
queueu(task_list_path, turn)



