import os
import re
import shutil
import subprocess
import sys
import tempfile


# Note sure why I did this. Test removal.
class Error(Exception):
    pass


def _write_script(cmds, script):
    """
    Write commands into a bash script
    """
    f = open(script, "w+")
    for cmd in cmds:
        print(cmd, file=f)
    f.close()


def decode(x):
    try:
        s = x.decode()
    except:
        return x
    return s


def cleanID(s):
    """
    Remove non alphanumeric characters from string.
    Replace whitespace with underscores.
    """
    s = re.sub(r"[^\w\s]", "", s)
    s = re.sub(r"\s+", "_", s)
    return s


def syscall(cmd, verbose=False):
    """
    Manage error handling when making syscalls
    """
    if verbose:
        print("Running command:", cmd, flush=True)
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        print(
            "The following command failed with exit code",
            error.returncode,
            file=sys.stderr,
        )
        print(cmd, file=sys.stderr)
        print("\nThe output was:\n", file=sys.stderr)
        print(error.output.decode(), file=sys.stderr)
        raise Error("Error running command:", cmd)
    if verbose:
        print(decode(output))


def run_cmd(cmds, verbose=False, tempDir=None, keeptemp=False):
    """
    Write and excute script
    """
    if not tempDir:
        tempDir = os.getcwd()
    tmpdir = tempfile.mkdtemp(prefix="tmp.", dir=tempDir)
    original_dir = os.getcwd()
    os.chdir(tmpdir)
    script = "run_jobs.sh"
    _write_script(cmds, script)
    syscall("bash " + script, verbose=verbose)
    os.chdir(original_dir)
    if not keeptemp:
        shutil.rmtree(tmpdir)
