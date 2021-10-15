#!/bin/bash
sudo apt-get update && 
	sudo apt-upgrade -y && 
	sudo apt-get install -y \
	build-essential \
	automake \
	libboost-dev-all \
	uuid-dev \
	libgpgme-dev \
	squashfs-tools \
	libseccomp-dev \
	wget \
	pkg-config \
	git \
	cryptsetup-bin

export VERSION=1.13.5 OS=linux ARCH=amd64 && \
    wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz && \
    sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz && \
    rm go$VERSION.$OS-$ARCH.tar.gz


echo 'export GOPATH=${HOME}/go' >> ~/.bashrc && \
    echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc && \
    source ~/.bashrc


export VERSION=3.8.0 && # adjust this as necessary \
    wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-ce-${VERSION}.tar.gz && \
    tar -xzf singularity-ce-${VERSION}.tar.gz && \
    cd singularity-ce-${VERSION}

export PATH=/usr/local/go/bin:$PATH

./mconfig && \
    make -C ./builddir && \
    sudo make -C ./builddir install

. /usr/local/etc/bash_completion.d/singularity


wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
bash ~/miniconda.sh -b -p $HOME/miniconda

mkdir -p $HOME/bin
cd $HOME/bin 
singularity build gromacs.simg docker:nvcr.io/hpc/gromacs:2021
echo '#!/bin/bash' > gmx
export SINGULARITYENV_DSSP=$DSSP
echo '$HOME/bin/gromacs.simg --nv gmx $*' >> gmx
chmod +x gmx
echo 'export PATH=$HOME/bin:$PATH' >> .bashrc	
