"""
Command execution wrappers with proper error handling.

Provides safe subprocess execution tools:
- Single and sequential command execution
- Temporary directory management

Emphasizes security by avoiding shell=True where possible
and using list-form commands.
"""

import logging
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import List, Optional, Sequence, Union


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
    cmds: Sequence[Union[str, List[str]]],
    verbose: bool = True,
    timeout: Optional[int] = None,
    cwd: Optional[Union[str, Path]] = None,
    stop_on_error: bool = True,
) -> List[Optional[subprocess.CompletedProcess]]:
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
    results: List[Optional[subprocess.CompletedProcess]] = []

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
) -> List[Optional[subprocess.CompletedProcess]]:
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
