function convert_dicoms_forSiemens(inpath)

dicom_dir = [inpath,filesep,'Dicoms'];
nifti_dir = [inpath, filesep, 'Niftis'];

% Go into the dicom dir 
cd(dicom_dir)

% Make a list of all the folders
folders = dir('0*');

% Loop through each folder name, and call dcm2niix
for f = 1:length(folders)
    unix(['dcm2niix -o ',nifti_dir,' -z y -f %f ./',folders(f).name],'-echo');
end

% Go back to the main subject directory
cd(inpath)

end