import subprocess
import os

# Assuming IGNITE_HOME is set in your environment variables
IGNITE_HOME = os.environ.get("IGNITE_HOME", "apache-ignite-2.15.0-bin")

def start_ignite_node(config_file):
    ignite_script = os.path.join(IGNITE_HOME, "bin", "ignite.sh")
    subprocess.Popen([ignite_script, config_file])

if __name__ == "__main__":
    # List of Ignite configuration files for each node
    config_files = [
        "apache-ignite-2.15.0-bin/config/ignite-config-1.xml",
        "apache-ignite-2.15.0-bin/config/ignite-config-2.xml",
        "apache-ignite-2.15.0-bin/config/ignite-config-3.xml",
        "apache-ignite-2.15.0-bin/config/ignite-config-4.xml",
    ]

    # Start each Ignite node in a separate process
    for config_file in config_files:
        start_ignite_node(config_file)

        
