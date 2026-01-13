"""
Command execution wrappers with proper error handling.

Provides safe subprocess execution tools:
- Single and sequential command execution
- Temporary directory management
- Script file generation and execution
- Legacy compatibility functions

Emphasizes security by avoiding shell=True where possible
and using list-form commands.
"""

import logging
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import List, Optional, Union


class CommandError(Exception):
    """
    Custom exception for command execution errors.

    Parameters
    ----------
    message : str
        Error message.
    cmd : str
        Command that failed.
    returncode : int
        Exit code from failed command.
    output : str, default ''
        Command output (stdout/stderr).
    """

    def __init__(self, message: str, cmd: str, returncode: int, output: str = ''):
        self.message = message
        self.cmd = cmd
        self.returncode = returncode
        self.output = output
        super().__init__(self.message)


def cleanID(s: str) -> str:
    """
    Remove non-alphanumeric characters and normalize whitespace in string.

    Parameters
    ----------
    s : str
        Input string to clean.

    Returns
    -------
    str
        Cleaned string with only alphanumeric characters and underscores.
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
    """
    Execute a system command with proper error handling and logging.

    Parameters
    ----------
    cmd : str or list of str
        Command to execute as string (requires shell=True) or list of arguments.
    verbose : bool, default False
        If True, logs command and output.
    timeout : int, optional
        Command timeout in seconds. None for no timeout.
    cwd : str or Path, optional
        Working directory for command execution.
    shell : bool, default False
        If True, executes command through shell (security risk with untrusted input).

    Returns
    -------
    subprocess.CompletedProcess
        Result object containing return code, stdout, and stderr.

    Raises
    ------
    CommandError
        If command fails (non-zero exit) or times out.

    Notes
    -----
    Prefer shell=False for security. Use list form of cmd when possible.
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
    """
    Execute multiple commands in sequence with error handling.

    Parameters
    ----------
    cmds : list of str or list of list of str
        List of commands to execute sequentially.
    verbose : bool, default True
        If True, logs progress and command output.
    timeout : int, optional
        Timeout in seconds applied to each command individually.
    cwd : str or Path, optional
        Working directory for all command executions.
    stop_on_error : bool, default True
        If True, stops execution on first failed command and raises error.
        If False, logs errors and continues with remaining commands.

    Returns
    -------
    list of subprocess.CompletedProcess or None
        Results from each command. Failed commands have None if stop_on_error=False.

    Raises
    ------
    CommandError
        If any command fails and stop_on_error is True.
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
    """
    Execute commands in a temporary directory with automatic cleanup.

    Parameters
    ----------
    cmds : list of str or list of list of str
        List of commands to execute in temporary directory.
    verbose : bool, default False
        If True, logs commands, output, and temp directory location.
    tempDir : str or Path, optional
        Parent directory for creating temporary directory. Uses cwd if None.
    keeptemp : bool, default False
        If True, preserves temporary directory after execution.
    timeout : int, optional
        Timeout in seconds applied to each command.

    Returns
    -------
    list of subprocess.CompletedProcess
        Results from all command executions.

    Raises
    ------
    CommandError
        If any command fails.

    Notes
    -----
    Modern replacement for legacy run_cmd() function.
    Always returns to original directory even if commands fail.
    Temporary directory named with 'tirmite_tmp_' prefix.
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
    """
    Write commands to an executable shell script file.

    Parameters
    ----------
    cmds : list of str
        List of command strings to write to script.
    script_path : str or Path
        Path for output script file.
    shell : str, default 'bash'
        Shell interpreter to use in shebang line.
    executable : bool, default True
        If True, sets file permissions to make script executable (chmod 755).

    Returns
    -------
    Path
        Path object for the created script file.

    Notes
    -----
    Adds shebang and error handling directives:
    - #!/bin/{shell}
    - set -euo pipefail (exit on error, undefined vars, pipe failures)
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
    """
    Execute a shell script file.

    Parameters
    ----------
    script_path : str or Path
        Path to script file to execute.
    verbose : bool, default False
        If True, logs command and output.
    timeout : int, optional
        Script timeout in seconds. None for no timeout.
    cwd : str or Path, optional
        Working directory for script execution.

    Returns
    -------
    subprocess.CompletedProcess
        Result object from script execution.

    Raises
    ------
    FileNotFoundError
        If script file doesn't exist.

    Notes
    -----
    Automatically sets execute permissions (chmod 755) before running.
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
    """
    DEPRECATED: Use run_command() instead.

    Parameters
    ----------
    cmd : str
        Command to execute.
    verbose : bool, default False
        Print verbose output.

    Returns
    -------
    None
        No return value. Legacy function for backwards compatibility.
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
    """
    DEPRECATED: Use run_cmd_in_tempdir() instead.

    Parameters
    ----------
    cmds : list of str
        Commands to execute.
    verbose : bool, default False
        Print verbose output.
    tempDir : str, optional
        Temporary directory path.
    keeptemp : bool, default False
        Keep temporary directory after execution.

    Returns
    -------
    None
        No return value. Legacy function for backwards compatibility.
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
    """
    DEPRECATED: Use write_script_file() instead.

    Parameters
    ----------
    cmds : list of str
        Commands to write to script.
    script : str
        Path to script file.

    Returns
    -------
    None
        No return value. Legacy function for backwards compatibility.
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
    """
    DEPRECATED: Modern subprocess uses text=True parameter.

    Parameters
    ----------
    x : bytes or str
        Value to decode.

    Returns
    -------
    str
        Decoded string.
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
    """Legacy exception class for backwards compatibility."""

    pass
