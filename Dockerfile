FROM python:3.11-slim

WORKDIR /sp_cont

RUN apt-get update && apt-get install -y build-essential git && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml ./
COPY src/ ./src/
COPY README.md ./

# Install the package using pip (setuptools will be used automatically)
RUN pip install .

# Entrypoint for the container -> If set all CLI commands can be directly used
ENTRYPOINT ["specimen"]
CMD ["-h"]