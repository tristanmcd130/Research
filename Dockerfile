FROM sagemath/sagemath
WORKDIR /usr/local/app

RUN sudo apt update && sudo apt install -y graphviz
RUN sage -pip install dot2tex