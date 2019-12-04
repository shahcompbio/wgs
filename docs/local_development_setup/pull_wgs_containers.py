import yaml
import sys
import subprocess

containers = yaml.load(open(sys.argv[1]), Loader=yaml.FullLoader)['wgs_containers']
cmd = "docker pull {}"
for name, container in containers.items():
    runcmd = cmd.format(container)
    process = subprocess.Popen(runcmd.split(), stdout = subprocess.PIPE)
    output, error = process.communicate()
    print ("pulling container {}".format(name), 
        "\noutput: ", output)
