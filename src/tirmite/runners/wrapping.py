import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import List, Optional, Union


class CommandError(Exception):
    """Custom exception for command execution errors."""

    def __init__(self, message: str, cmd: str, returncode: int, output: str = ''):
        self.message = message
        self.cmd = cmd
        self.returncode = returncode
        self.output = output
        super().__init__(self.message)


def cleanID(s: str) -> str:
    """Remove non alphanumeric characters from string.
    Replace whitespace with underscores.

    Args:
        s: Input string to clean

    Returns:
        str: Cleaned string with only alphanumeric chars and underscores

    """
    s = re.sub(r'[^\w\s]', '', s)
    s = re.sub(r'\s+', '_', s)
    return s


def run_command(
    cmd: Union[str, List[str]],
    verbose: bool = False,
    timeout: Optional[int] = None,
    cwd: Optional[Union[str, Path]] = None,
    shell: bool = False,
) -> subprocess.CompletedProcess:
    """Execute a command using subprocess with proper error handling.

    Args:
        cmd: Command to execute (string or list of arguments)
        verbose: Print command and output if True
        timeout: Command timeout in seconds (None for no timeout)
        cwd: Working directory for command execution
        shell: Whether to use shell execution (discouraged for security)

    Returns:
        subprocess.CompletedProcess: Result of command execution

    Raises:
        CommandError: If command fails or times out

    """
    if verbose:
        cmd_str = cmd if isinstance(cmd, str) else ' '.join(cmd)
        logging.info(f'Running command: {cmd_str}')

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=cwd,
            shell=shell,
            check=False,  # Don't raise CalledProcessError, handle manually
        )

        if result.returncode != 0:
            cmd_str = cmd if isinstance(cmd, str) else ' '.join(cmd)
            error_msg = f'Command failed with exit code {result.returncode}'

            # Combine stdout and stderr for error output
            output = ''
            if result.stdout:
                output += f'STDOUT:\n{result.stdout}\n'
            if result.stderr:
                output += f'STDERR:\n{result.stderr}\n'

            raise CommandError(
                message=error_msg,
                cmd=cmd_str,
                returncode=result.returncode,
                output=output,
            )

        if verbose and result.stdout:
            logging.info(f'Command output:\n{result.stdout}')

        return result

    except subprocess.TimeoutExpired as err:
        cmd_str = cmd if isinstance(cmd, str) else ' '.join(cmd)
        raise CommandError(
            message=f'Command timed out after {timeout} seconds',
            cmd=cmd_str,
            returncode=-1,
        ) from err
    except Exception as e:
        cmd_str = cmd if isinstance(cmd, str) else ' '.join(cmd)
        raise CommandError(
            message=f'Error executing command: {str(e)}', cmd=cmd_str, returncode=-1
        ) from e


def run_commands_sequential(
    cmds: List[Union[str, List[str]]],
    verbose: bool = True,
    timeout: Optional[int] = None,
    cwd: Optional[Union[str, Path]] = None,
    stop_on_error: bool = True,
) -> List[subprocess.CompletedProcess]:
    """Execute multiple commands sequentially.

    Args:
        cmds: List of commands to execute
        verbose: Print commands and output if True
        timeout: Timeout per command in seconds
        cwd: Working directory for command execution
        stop_on_error: Stop execution if a command fails

    Returns:
        List of subprocess.CompletedProcess results

    Raises:
        CommandError: If any command fails and stop_on_error is True

    """
    results = []

    for i, cmd in enumerate(cmds):
        try:
            logging.info(f'Executing command {i + 1}/{len(cmds)}')

            result = run_command(
                cmd=cmd,
                verbose=verbose,
                timeout=timeout,
                cwd=cwd,
                shell=isinstance(cmd, str),  # Use shell only for string commands
            )
            results.append(result)

        except CommandError as e:
            if stop_on_error:
                logging.error(f'Command {i + 1} failed: {e.message}')
                raise
            else:
                logging.warning(f'Command {i + 1} failed (continuing): {e.message}')
                results.append(None)  # Placeholder for failed command
                continue

    return results


def run_cmd_in_tempdir(
    cmds: List[Union[str, List[str]]],
    verbose: bool = False,
    tempDir: Optional[Union[str, Path]] = None,
    keeptemp: bool = False,
    timeout: Optional[int] = None,
) -> List[subprocess.CompletedProcess]:
    """Execute commands in a temporary directory with automatic cleanup.

    This is the modern replacement for the original run_cmd function.

    Args:
        cmds: List of commands to execute
        verbose: Print commands and output if True
        tempDir: Parent directory for temporary directory (default: current dir)
        keeptemp: Keep temporary directory after execution
        timeout: Timeout per command in seconds

    Returns:
        List of subprocess.CompletedProcess results

    Raises:
        CommandError: If any command fails

    """
    if tempDir is None:
        tempDir = os.getcwd()

    tempDir = Path(tempDir)
    original_dir = Path.cwd()

    # Create temporary directory
    tmpdir = None
    try:
        tmpdir = tempfile.mkdtemp(prefix='tirmite_tmp_', dir=tempDir)
        tmpdir_path = Path(tmpdir)

        if verbose:
            logging.info(f'Working in temporary directory: {tmpdir_path}')

        # Execute commands in temporary directory
        results = run_commands_sequential(
            cmds=cmds,
            verbose=verbose,
            timeout=timeout,
            cwd=tmpdir_path,
            stop_on_error=True,
        )

        return results

    finally:
        # Always return to original directory
        os.chdir(original_dir)

        # Clean up temporary directory if requested
        if tmpdir and not keeptemp:
            try:
                shutil.rmtree(tmpdir)
                if verbose:
                    logging.info(f'Cleaned up temporary directory: {tmpdir}')
            except OSError as e:
                logging.warning(f'Failed to remove temporary directory {tmpdir}: {e}')


def write_script_file(
    cmds: List[str],
    script_path: Union[str, Path],
    shell: str = 'bash',
    executable: bool = True,
) -> Path:
    """Write commands to a script file with proper headers and permissions.

    Args:
        cmds: List of commands to write
        script_path: Path for the script file
        shell: Shell interpreter (default: 'bash')
        executable: Make script executable

    Returns:
        Path: Path to created script file

    """
    script_path = Path(script_path)

    # Ensure parent directory exists
    script_path.parent.mkdir(parents=True, exist_ok=True)

    # Write script with proper shebang
    with open(script_path, 'w') as f:
        f.write(f'#!/bin/{shell}\n')
        f.write('set -euo pipefail\n')  # Exit on error, undefined vars, pipe failures
        f.write('\n')

        for cmd in cmds:
            f.write(f'{cmd}\n')

    # Make executable if requested
    if executable:
        script_path.chmod(0o755)

    return script_path


def run_script_file(
    script_path: Union[str, Path],
    verbose: bool = False,
    timeout: Optional[int] = None,
    cwd: Optional[Union[str, Path]] = None,
) -> subprocess.CompletedProcess:
    """Execute a script file.

    Args:
        script_path: Path to script file
        verbose: Print command and output if True
        timeout: Script timeout in seconds
        cwd: Working directory for script execution

    Returns:
        subprocess.CompletedProcess: Result of script execution

    """
    script_path = Path(script_path)

    if not script_path.exists():
        raise FileNotFoundError(f'Script file not found: {script_path}')

    # Make sure script is executable
    script_path.chmod(0o755)

    return run_command(
        cmd=[str(script_path)], verbose=verbose, timeout=timeout, cwd=cwd, shell=False
    )


# Legacy functions for backwards compatibility
def syscall(cmd: str, verbose: bool = False) -> None:
    """DEPRECATED: Use run_command() instead.
    Legacy function maintained for backwards compatibility.
    """
    import warnings

    warnings.warn(
        'syscall() is deprecated. Use run_command() instead.',
        DeprecationWarning,
        stacklevel=2,
    )

    try:
        run_command(cmd, verbose=verbose, shell=True)
    except CommandError as e:
        # Convert to original Error for backwards compatibility
        raise Error(f'Error running command: {e.cmd}') from e


def run_cmd(
    cmds: List[str],
    verbose: bool = False,
    tempDir: Optional[str] = None,
    keeptemp: bool = False,
) -> None:
    """DEPRECATED: Use run_cmd_in_tempdir() instead.
    Legacy function maintained for backwards compatibility.
    """
    import warnings

    warnings.warn(
        'run_cmd() is deprecated. Use run_cmd_in_tempdir() instead.',
        DeprecationWarning,
        stacklevel=2,
    )

    try:
        run_cmd_in_tempdir(cmds, verbose, tempDir, keeptemp)
    except CommandError as e:
        # Convert to original Error for backwards compatibility
        raise Error(f'Error running commands: {e.message}') from e


def _write_script(cmds: List[str], script: str) -> None:
    """DEPRECATED: Use write_script_file() instead.
    Legacy function maintained for backwards compatibility.
    """
    import warnings

    warnings.warn(
        '_write_script() is deprecated. Use write_script_file() instead.',
        DeprecationWarning,
        stacklevel=2,
    )

    with open(script, 'w') as f:
        for cmd in cmds:
            print(cmd, file=f)


def decode(x):
    """DEPRECATED: Modern subprocess uses text=True parameter.
    Legacy function maintained for backwards compatibility.
    """
    import warnings

    warnings.warn(
        'decode() is deprecated. Use text=True in subprocess.run() instead.',
        DeprecationWarning,
        stacklevel=2,
    )

    try:
        return x.decode()
    except (AttributeError, UnicodeDecodeError):
        return x


# Keep original Error class for backwards compatibility
class Error(Exception):
    pass
