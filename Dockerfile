FROM centos:centos7.3.1611
MAINTAINER lucioric2000@hotmail.com

ADD install_qiaseq_dna.bash /srv
RUN yum -y install sudo
RUN chmod +x /srv/install_qiaseq_dna.bash
RUN cd /srv; ./install_qiaseq_dna.bash

EXPOSE 80

CMD source /srv/conda/bin/conda activate base