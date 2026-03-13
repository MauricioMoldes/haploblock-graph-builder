FROM python:3.11-slim

# metadata
LABEL maintainer="Mauricio"
LABEL description="Haploblock graph builder for synthetic node expansion pipeline"

# working directory
WORKDIR /app

# copy project
COPY haploblocks_to_graph.py /app/
COPY requirements.txt /app/

# creates dirs 
RUN mkdir /data /results

# install dependencies
RUN pip install --no-cache-dir -r requirements.txt

# default command
ENTRYPOINT ["python", "/app/haploblocks_to_graph.py"]
