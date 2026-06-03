#!/bin/bash
#This script installs MEME (including dependencies)

##
##Make software directories
mkdir -p ../../software/cpan
mkdir -p ../../software/meme

##
##Download meme
wget https://meme-suite.org/meme/meme-software/5.5.5/meme-5.5.5.tar.gz -P ../../software/meme/
tar zxf ../../software/meme/meme-5.5.5.tar.gz -C ../../software/meme/

##
##Install perl libraries
cpanm -l ../../software/cpan/ -L ../../software/cpan/ File::Which HTML::Template JSON XML::Simple Sys::Info Math::CDF

##
##Install MEME
cd ../../software/meme/meme-5.5.5/
./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt
make
make test
make install

##
##Add path
echo 'export PATH=$HOME/meme/bin:$PATH' >> ~/.bashrc
