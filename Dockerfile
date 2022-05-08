# The first line to add
# syntax=docker/dockerfile:1

FROM python:3.8-slim-buster
WORKDIR /app
# Copy the requirements.txt file into the working directory
COPY requirements.txt requirements.txt
# Install python pip requirements 
RUN pip3 install -r requirements.txt
# Copy all the files 
COPY . .

CMD []


