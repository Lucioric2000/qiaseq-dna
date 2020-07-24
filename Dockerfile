FROM centos:centos7.3.1611
MAINTAINER lucioric2000@hotmail.com
ENV version 15.5
ENV installerloc=/root
ENV installer=install_qiaseq-dna-v${version}.bash
ENV tarfile=qiaseq-dna-${version}.tar.gz
VOLUME /srv/qgen/data
RUN ls /srv/qgen/data
RUN mkdir -p ${installerloc}
RUN yum -y update && yum -y install sudo
ADD ${installer} ${installerloc}/
ADD $tarfile ${installerloc}/
RUN chmod +x ${installerloc}/${installer}
RUN cd ${installerloc} && ./${installer}


