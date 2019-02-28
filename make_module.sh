git clone https://github.com/lukaselflein/smamp

mkdir ~/modulefiles

cd smamp
echo "#%Module1.0" > ~/modulefiles/smamp
echo "prepend-path PYTHONPATH $(pwd)" >> ~/modulefiles/smamp

module use ~/modulefiles
module load smamp
