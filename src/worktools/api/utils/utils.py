import subprocess


def run_shell_command(command, output_file_path=None):
    """Execute a shell command and store the output to a file.

    Args:
        command (str): The shell command to execute.
        output_file_path (str): Path to the file where the command output will be stored.
    """
    output = ""
    try:
        # Execute the command
        output = subprocess.run(
            command, check=True, shell=True, capture_output=True
        ).stdout.decode("utf-8")
        # Write the output to the specified file
        output_file_path = "/home/rg972/project/intersect_commmand.txt"
        with open(output_file_path, "w") as file:
            file.write(output)
    except subprocess.CalledProcessError as e:
        output = f"Error executing command: {e}"
        raise Exception(output)

    return output


def bgzip_file(file_path, output_file_path=None):
    if output_file_path is None:
        output_file_path = file_path + ".gz"

    bgzip_cmd = f"bgzip -c {file_path} > {output_file_path}"
    run_shell_command(bgzip_cmd)

    return output_file_path


def is_bgzipped(file_path):
    with open(file_path, "rb") as file:
        # Read the gzip signature from the start of the file
        gzip_signature = file.read(2)
        if gzip_signature != b"\x1f\x8b":
            return False

        # Skip the next 8 bytes to move to the XFL (extra flags) byte
        file.seek(8, 1)

        # Read the extra field length (XLEN), which is 2 bytes, little-endian
        xlen_bytes = file.read(2)
        if len(xlen_bytes) < 2:
            return False

        # Unpack XLEN to get the length of the extra field
        xlen = int.from_bytes(xlen_bytes, byteorder="little")

        # Read the extra field itself
        extra_field = file.read(xlen)

        # BGZF uses a specific subfield in the extra field with SI1 = 66 (0x42) and SI2 = 67 (0x43)
        # indicating a BGZF block if present
        if b"\x42\x43" in extra_field:
            return True

    return False


def get_num_lines_in_file(file_path):
    with open(file_path, "r") as file:
        return sum(1 for line in file)
