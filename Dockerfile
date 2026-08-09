FROM python:3.11-slim

WORKDIR /app

# FPSim2 / RDKit / pyarrow wheels are self-contained; no system chem libs needed.
COPY pyproject.toml README.md ./
COPY build_index.py server.py ./
RUN pip install --no-cache-dir .

COPY entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

EXPOSE 8080
# Usage:
#   docker run ... build --input <parquet|s3 glob> --out /index
#   docker run -p 8080:8080 ... serve --index /index [--mode in-memory|on-disk]
ENTRYPOINT ["entrypoint.sh"]
CMD ["serve", "--index", "/index"]
