"""
Environment manager for BaseBuddy GUI
Handles automatic switching between ARM and x86 environments
"""

import subprocess
import sys
import os
import logging
from pathlib import Path
from typing import Optional, Dict, Any

logger = logging.getLogger(__name__)

class EnvironmentManager:
    """Manages conda environment switching for GUI operations"""
    
    def __init__(self):
        self.current_env = os.environ.get('CONDA_DEFAULT_ENV', '')
        self.conda_base = self._get_conda_base()
        self.environments = {
            'arm': 'basebuddy-arm',
            'x86': 'basebuddy-x86'
        }
        
    def _get_conda_base(self) -> str:
        """Get conda base directory"""
        try:
            result = subprocess.run(['conda', 'info', '--base'], 
                                  capture_output=True, text=True, check=True)
            return result.stdout.strip()
        except:
            return os.path.expanduser('~/miniforge3')
    
    def check_tool_availability(self, tool: str, env: str = None) -> bool:
        """Check if a tool is available in specified environment"""
        if env:
            cmd = ['conda', 'run', '-n', self.environments.get(env, env), 'which', tool]
        else:
            cmd = ['which', tool]
            
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            return result.returncode == 0
        except:
            return False
    
    def run_in_environment(self, cmd: list, env: str = 'arm') -> subprocess.CompletedProcess:
        """Run a command in specified conda environment"""
        env_name = self.environments.get(env, env)
        
        # Use conda run to execute in specific environment
        full_cmd = ['conda', 'run', '-n', env_name] + cmd
        
        logger.info(f"Running in {env_name}: {' '.join(cmd)}")
        
        return subprocess.run(full_cmd, capture_output=True, text=True)
    
    def get_tool_environment(self, tool: str) -> Optional[str]:
        """Determine which environment a tool should run in"""
        # Tools that require x86 environment
        x86_tools = {
            'addsnv.py', 'addindel.py', 'addsv.py',
            'exonerate', 'velvet', 'velveth', 'velvetg',
            'wgsim', 'nanosim', 'nanosim-h'
        }
        
        if tool in x86_tools:
            return 'x86'
        return 'arm'
    
    def validate_environments(self) -> Dict[str, bool]:
        """Check which environments are properly installed"""
        status = {}
        for name, env in self.environments.items():
            try:
                cmd = ['conda', 'env', 'list']
                result = subprocess.run(cmd, capture_output=True, text=True)
                status[name] = env in result.stdout
            except:
                status[name] = False
        return status


# Singleton instance
env_manager = EnvironmentManager()


def run_with_auto_env(command: list, **kwargs) -> subprocess.CompletedProcess:
    """
    Run a command with automatic environment selection
    
    Args:
        command: Command as list of strings
        **kwargs: Additional arguments for subprocess.run
        
    Returns:
        CompletedProcess instance
    """
    if not command:
        raise ValueError("Command cannot be empty")
    
    # Determine which tool is being called
    tool = os.path.basename(command[0])
    
    # Get appropriate environment
    env = env_manager.get_tool_environment(tool)
    
    # Run in selected environment
    if env and env != 'arm':  # Only use conda run for x86 env
        full_cmd = ['conda', 'run', '-n', env_manager.environments[env]] + command
        logger.info(f"Running {tool} in {env} environment")
    else:
        full_cmd = command
        
    return subprocess.run(full_cmd, **kwargs)