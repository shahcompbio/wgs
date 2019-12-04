import yaml
import sys
import subprocess
import argparse

def pull_containers(containers):
    '''
    pulls all docker containers from yaml
    file config.

    :param containers: yaml with docker containers
    '''
    cmd = "docker pull {}"
    for name, container in containers.items():
        runcmd = cmd.format(container)
        process = subprocess.Popen(runcmd.split(), stdout = subprocess.PIPE)
        output, error = process.communicate()
        print ("pulling container {}".format(name), 
            "\noutput: ", output)

def main(): 

    parser = argparse.ArgumentParser()
    parser.add_argument('containers_yaml', 
        help='yaml containing the wgs docker containers. See containers.yaml')
    
    args = parser.parse_args()
    
    containers_yaml = vars(args)["containers_yaml"]
    containers = yaml.load(open(containers_yaml), Loader=yaml.FullLoader)['wgs_containers']
    
    pull_containers(containers)

if __name__ == "__main__":
    main()

