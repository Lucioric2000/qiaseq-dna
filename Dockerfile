FROM centos:centos7.3.1611
MAINTAINER lucioric2000@hotmail.com
ENV version=14.1
ENV installer=install_qiaseq-dna-v${version}.bash
ENV tarfileqiaseq-dna-v${version}.bash

RUN mkdir -p /srv/qgen
ADD ${installer}.bash /srv/qgen
RUN yum -y install sudo
RUN chmod +x /srv/qgen/${installer}
RUN cd /srv/qgen; ./${installer}

EXPOSE 80

CMD source /root/bin/conda activate base