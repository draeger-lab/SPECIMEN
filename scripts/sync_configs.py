import os

# Define paths relative to the root directory
CONFIG_DIR = os.path.join("src", "specimen", "data", "config")
DEFAULT_CONFIG = os.path.join(CONFIG_DIR, "hqtb_config_default.yaml")
ADVANCED_CONFIG = os.path.join(CONFIG_DIR, "hqtb_advanced_config_expl.yaml")
BASIC_CONFIG = os.path.join(CONFIG_DIR, "hqtb_basic_config_expl.yaml")

# Developer tags to strip out for non-default configs
DEV_TAGS = ["# @IDEA", "# @TODO", "# @DEV"]

# Blocks to entirely skip for the basic "quick-and-dirty" configuration
ADVANCED_BLOCKS = [
    "refinement_cleanup:",
    "refinement_smoothing:",
    "GeneGapFiller:",
    "GeneGapFiller parameters:",
    "media_gap:",
    "mcc:",
    "egc:",
]

def clean_dev_tags(line: str) -> str:
    """Strips dev tags from a line. Returns the cleaned line, or empty string if it was standalone."""
    for tag in DEV_TAGS:
        if tag in line:
            # Split at the tag and keep the left side
            cleaned = line.split(tag)[0].rstrip()
            # If there's content left (like a YAML key), add the newline back
            if cleaned:
                return cleaned + "\n"
            return ""
    return line

def generate_advanced(lines: list[str]) -> list[str]:
    """Generates advanced config by stripping dev notes."""
    advanced_lines = []
    for line in lines:
        cleaned = clean_dev_tags(line)
        if cleaned:
            advanced_lines.append(cleaned)
    return advanced_lines

def generate_basic(lines: list[str]) -> list[str]:
    """Generates basic config by stripping dev notes and advanced blocks."""
    basic_lines = []
    skip_block = False
    current_indent = 0

    for line in lines:
        cleaned = clean_dev_tags(line)
        if not cleaned:
            continue

        stripped = cleaned.lstrip()
        indent = len(cleaned) - len(stripped)

        # Handle skipping nested blocks
        if skip_block:
            if stripped and indent <= current_indent:
                skip_block = False
            else:
                continue

        # Check if we should start skipping a new block
        if any(stripped.startswith(key) for key in ADVANCED_BLOCKS):
            skip_block = True
            current_indent = indent
            continue

        basic_lines.append(cleaned)

    return basic_lines


def main():
    print("Syncing SPECIMEN YAML configurations...")
    if not os.path.exists(DEFAULT_CONFIG):
        print(
            f"Error: Could not find {DEFAULT_CONFIG}. Run this script from the project root."
        )
        return

    with open(DEFAULT_CONFIG, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # Generate and write Advanced Config
    advanced_lines = generate_advanced(lines)
    with open(ADVANCED_CONFIG, "w", encoding="utf-8") as f:
        f.writelines(advanced_lines)
    print(f"  [+] Successfully generated {ADVANCED_CONFIG}")

    # Generate and write Basic Config
    basic_lines = generate_basic(lines)
    with open(BASIC_CONFIG, "w", encoding="utf-8") as f:
        f.writelines(basic_lines)
    print(f"  [+] Successfully generated {BASIC_CONFIG}")
    print("Done!")


if __name__ == "__main__":
    main()
