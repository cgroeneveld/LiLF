find . -type f -name '*.py' -exec sed -i 's/np\.float16/float/g' {} \;
find . -type f -name '*.py' -exec sed -i 's/np\.float/float/g' {} \;
find . -type f -name '*.py' -exec sed -i 's/np\.int/int/g' {} \;
find . -type f -name '*.py' -exec sed -i 's/np\.complex/complex/g' {} \;
find . -type f -name '*.py' -exec sed -i 's/np\.bool/bool/g' {} \;
