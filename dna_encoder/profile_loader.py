"""
Profile loader for DNA generator validation profiles.

Simple profile loading with two-level priority:
1. default_profiles.json (shipped with package)
2. user_profiles.json (in current working directory, optional)

User profiles override defaults with same name.
"""

import json
import logging
from pathlib import Path
from typing import Dict, Any, Optional

logger = logging.getLogger(__name__)


class ProfileLoader:
    """Loader for validation profiles from JSON files."""

    def __init__(self):
        """Initialize profile loader."""
        # Default profiles - shipped with package
        self.package_dir = Path(__file__).parent
        self.default_profiles_path = self.package_dir / 'default_profiles.json'

        # User profiles - in current working directory
        self.user_profiles_path = Path.cwd() / 'user_profiles.json'

        self._cache: Optional[Dict[str, Any]] = None

    def load_profiles(self, force_reload: bool = False) -> Dict[str, Any]:
        """
        Load all validation profiles.

        Priority:
        1. Load default_profiles.json (required)
        2. Load user_profiles.json (optional, overrides defaults)

        Args:
            force_reload: Force reload from disk even if cached

        Returns:
            Dictionary mapping profile names to profile configurations

        Raises:
            FileNotFoundError: If default_profiles.json not found
            ValueError: If default_profiles.json is invalid
        """
        if self._cache and not force_reload:
            return self._cache

        # 1. Load defaults (REQUIRED)
        if not self.default_profiles_path.exists():
            raise FileNotFoundError(
                f"Default profiles not found: {self.default_profiles_path}\n"
                f"This file should be shipped with the package."
            )

        default_data = self._load_json_file(self.default_profiles_path)
        if not default_data or 'profiles' not in default_data:
            raise ValueError(
                f"Invalid default_profiles.json: missing 'profiles' key"
            )

        profiles = default_data['profiles'].copy()
        logger.debug(f"Loaded {len(profiles)} default profiles")

        # 2. Load user profiles (OPTIONAL)
        if self.user_profiles_path.exists():
            user_data = self._load_json_file(self.user_profiles_path)
            if user_data and 'profiles' in user_data:
                user_profiles = user_data['profiles']
                # Merge: user profiles override defaults
                for name, profile in user_profiles.items():
                    if name in profiles:
                        logger.info(f"User profile '{name}' overrides default")
                    profiles[name] = profile
                logger.info(f"Loaded {len(user_profiles)} user profiles from {self.user_profiles_path}")
            else:
                logger.warning(f"Invalid user_profiles.json - ignoring")
        else:
            logger.debug(f"No user profiles found at {self.user_profiles_path}")

        # 3. Validate all profiles
        validated_profiles = {}
        for name, profile in profiles.items():
            if self._validate_profile(name, profile):
                validated_profiles[name] = profile
            else:
                logger.warning(f"Skipping invalid profile: {name}")

        if not validated_profiles:
            raise ValueError("No valid profiles found!")

        self._cache = validated_profiles
        return validated_profiles

    def get_profile(self, name: str) -> Optional[Dict[str, Any]]:
        """
        Get a specific profile by name.

        Args:
            name: Profile name

        Returns:
            Profile configuration or None if not found
        """
        profiles = self.load_profiles()
        return profiles.get(name)

    def list_profiles(self) -> Dict[str, str]:
        """
        List all available profiles with descriptions.

        Returns:
            Dictionary mapping profile names to descriptions
        """
        profiles = self.load_profiles()
        return {
            name: profile.get('description', 'No description')
            for name, profile in profiles.items()
        }

    def _load_json_file(self, path: Path) -> Optional[Dict[str, Any]]:
        """
        Load JSON file with error handling.

        Args:
            path: Path to JSON file

        Returns:
            Parsed JSON data or None if error
        """
        try:
            with open(path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            return data

        except json.JSONDecodeError as e:
            logger.error(f"Invalid JSON in {path}: {e}")
            return None
        except Exception as e:
            logger.error(f"Error loading {path}: {e}")
            return None

    def _validate_profile(self, name: str, profile: Dict[str, Any]) -> bool:
        """
        Validate profile structure.

        Args:
            name: Profile name
            profile: Profile configuration

        Returns:
            True if valid, False otherwise
        """
        # Required keys
        if 'rules' not in profile:
            logger.error(f"Profile '{name}' missing 'rules' section")
            return False

        if 'params' not in profile:
            logger.error(f"Profile '{name}' missing 'params' section")
            return False

        # Validate rules (must be dict of booleans)
        rules = profile['rules']
        if not isinstance(rules, dict):
            logger.error(f"Profile '{name}': 'rules' must be a dictionary")
            return False

        valid_rule_names = {
            'gc_content', 'melting_temperature', 'hairpin_structures',
            'homodimer_structures', 'homopolymer_runs',
            'dinucleotide_repeats', 'three_prime_stability'
        }

        for rule_name, value in rules.items():
            if rule_name not in valid_rule_names:
                logger.warning(f"Profile '{name}': unknown rule '{rule_name}'")
            if not isinstance(value, bool):
                logger.error(f"Profile '{name}': rule '{rule_name}' must be boolean")
                return False

        # Validate params (must be dict of numbers)
        params = profile['params']
        if not isinstance(params, dict):
            logger.error(f"Profile '{name}': 'params' must be a dictionary")
            return False

        valid_param_names = {
            'min_gc', 'max_gc', 'min_tm', 'max_tm',
            'max_hairpin_tm', 'max_homodimer_tm',
            'max_homopolymer_length', 'max_dinucleotide_repeats',
            'max_3prime_gc'
        }

        for param_name, value in params.items():
            if param_name not in valid_param_names:
                logger.warning(f"Profile '{name}': unknown param '{param_name}'")
            if not isinstance(value, (int, float)):
                logger.error(f"Profile '{name}': param '{param_name}' must be numeric")
                return False

        return True


# Global singleton instance
_loader: Optional[ProfileLoader] = None


def get_profile_loader() -> ProfileLoader:
    """Get global ProfileLoader instance (singleton)."""
    global _loader
    if _loader is None:
        _loader = ProfileLoader()
    return _loader


def load_profile(name: str) -> Optional[Dict[str, Any]]:
    """
    Convenience function to load a profile by name.

    Args:
        name: Profile name

    Returns:
        Profile configuration or None if not found
    """
    return get_profile_loader().get_profile(name)


def list_available_profiles() -> Dict[str, str]:
    """
    Convenience function to list all available profiles.

    Returns:
        Dictionary mapping profile names to descriptions
    """
    return get_profile_loader().list_profiles()
