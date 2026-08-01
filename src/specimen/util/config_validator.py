import os
import logging
from typing import Any, Dict, List, Optional, Union

logger = logging.getLogger(__name__)


class ConfigValidationError(Exception):
    """Custom exception for configuration validation errors."""

    pass


class SchemaField:
    """Defines the expected structure, type, and constraints of a config field."""

    def __init__(
        self,
        expected_type: type | tuple,
        default: Any = None,
        required: bool = False,
        min_val: Optional[Union[int, float]] = None,
        max_val: Optional[Union[int, float]] = None,
        is_path: bool = False,
        path_must_exist: bool = False,
        allowed_values: Optional[List[Any]] = None,
    ):
        self.expected_type = expected_type
        self.default = default
        self.required = required
        self.min_val = min_val
        self.max_val = max_val
        self.is_path = is_path
        self.path_must_exist = path_must_exist
        self.allowed_values = allowed_values


def validate_config(
    config: Dict[str, Any], schema: Dict[str, Any], path: str = ""
) -> Dict[str, Any]:
    """
    Recursively validates a configuration dictionary against a defined schema.
    Replaces verbose YAML traversal and USER/_USER_ string checks.
    """
    validated = {}

    for key, rule in schema.items():
        current_path = f"{path}.{key}" if path else key
        value = config.get(key, rule.default if isinstance(rule, SchemaField) else None)

        # Handle nested dictionaries
        if isinstance(rule, dict):
            if not isinstance(value, dict):
                if value is None:
                    value = {}
                else:
                    raise ConfigValidationError(
                        f"Expected a dictionary at '{current_path}', got {type(value).__name__}."
                    )
            validated[key] = validate_config(value, rule, current_path)
            continue

        if not isinstance(rule, SchemaField):
            raise ConfigValidationError(
                f"Invalid schema definition at '{current_path}'."
            )

        # 1. Check Required fields (Replaces _USER_ and USER logic)
        if value is None or value == "_USER_" or value == "USER":
            if rule.required:
                raise ConfigValidationError(
                    f"Missing required configuration parameter: '{current_path}'"
                )

            if value == "_USER_":
                logger.error(
                    f"Deprecated '_USER_' found at '{current_path}'. Please provide a value or remove it."
                )
                raise ConfigValidationError(
                    f"Missing required parameter: '{current_path}'"
                )

            if value == "USER":
                logger.warning(
                    f"Deprecated 'USER' keyword found at '{current_path}'. Using default: {rule.default}"
                )

            validated[key] = rule.default
            continue

        # Clean up legacy 'USER' keyword
        if value == "USER":
            logger.warning(
                f"Deprecated 'USER' keyword found at '{current_path}'. Using default: {rule.default}"
            )
            value = rule.default

        # 2. Type Checking
        if value is not None and not isinstance(value, rule.expected_type):
            raise ConfigValidationError(
                f"Type mismatch at '{current_path}': Expected {rule.expected_type}, got {type(value).__name__}."
            )

        # 3. Numeric Bounds Checking
        if isinstance(value, (int, float)):
            if rule.min_val is not None and value < rule.min_val:
                raise ConfigValidationError(
                    f"Value at '{current_path}' ({value}) is below minimum allowed ({rule.min_val})."
                )
            if rule.max_val is not None and value > rule.max_val:
                raise ConfigValidationError(
                    f"Value at '{current_path}' ({value}) is above maximum allowed ({rule.max_val})."
                )

        # 4. Allowed Values Checking (Strings, etc.)
        if rule.allowed_values is not None and value not in rule.allowed_values:
            raise ConfigValidationError(
                f"Value at '{current_path}' ({value}) must be one of {rule.allowed_values}."
            )

        # 5. Path Validation
        if rule.is_path and value is not None:
            value = str(value)
            if rule.path_must_exist and not os.path.exists(value):
                raise ConfigValidationError(
                    f"Path specified at '{current_path}' does not exist: {value}"
                )

        validated[key] = value

    return validated
