FROM ubuntu:latest

RUN  apt-get update \
  && apt-get install -y wget curl
  

RUN cd /usr/local/bin && wget https://github.com/Illumina/ExpansionHunter/releases/download/v5.0.0/ExpansionHunter-v5.0.0-linux_x86_64.tar.gz \
    && tar -xvf ExpansionHunter-v5.0.0-linux_x86_64.tar.gz \
    && mv /usr/local/bin/ExpansionHunter-v5.0.0-linux_x86_64/bin/ExpansionHunter /usr/local/bin \
    && rm -f /usr/local/bin/ExpansionHunter-v5.0.0-linux_x86_64.tar.gz \
    && rm -rf /usr/local/bin/ExpansionHunter-v5.0.0-linux_x86_64/example