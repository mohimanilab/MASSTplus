# /projects/mohimanilab/OpenMS-2.3.0/bin/FileConverter -in 000001527_RG8_01_6306.mzML -out 000001527_RG8_01_6306.mgf

source /projects/mohimanilab/.bashrc;

echo 'Converting all Spectras from .mzML to .mgf';

for f in data/*.mzML; do
    # do some stuff here with "$f"
    
    echo "=============================================";
    echo "Converting ${f}...";
    /projects/mohimanilab/OpenMS-2.3.0/bin/FileConverter -in "${f}" -out "${f}.mgf"
    echo "Done converting ${f}...";
done

echo 'Finished converting all Spectras from .mzML to .mgf!';
