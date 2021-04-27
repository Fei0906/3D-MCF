%================================================================
%   This file is used to set parameters, input files, and output
%   files for 3D_MCF function. Please first copy this file to your
%   working directory before editing. Then modify the value of the 
%   parameters as you need. However, please DO NOT modify the 
%   names of parameters.
%
%   Fei Liu, Oct 2019
%================================================================


%================================================================
%	!SAVE IMPROTANT TEMPORARY DATA BEFORE RUN 3D_MCF!
%================================================================
clear;


%================================================================
%   Input and output file. To unify the interface, all the input 
%   and output data are formed in MATLAB data structure (stored
%   in *.mat file). With default value, you should provide a file
%   called 'inputData.mat' which contains the variables of 'wPh', 
%   'wgt', 'ifgs', and 'coor'. And the data is stored in these
%   variables.
%================================================================

%   The input file name, which contains all the required data in it.
INPUTFILE='input_3DMCF_data_3.mat';

%   The flag which indicates the input data type, where 'C' and 
%   'D' stand for the continuous data and the discrete data
%   respectively.
DATATYPE='D';

%   Wrapped phases, which are stored in INPUTFILE. When data type
%   is 'C', it should be a size of M*N*P matrix, where M, N, and 
%   P stand for the row, the column, and the number of the 
%   interferograms respectively. And when data type is 'D', it
%   should be a size of S*P, where S stands for the number of 
%   pixels in individual interferogram.
WRAPPEDPHASES='wPh';

%   Weight coefficients, which are stored in INPUTFILE. It should
%   be the same size of WRAPPEDPHASES, and stands for the weight 
%   value used in our algorithm.
WEIGHTCOEFFICIENTS='wgt';

%   The connections of SAR images, which are stored in INPUTFILE.
%   It should be a size of P*2 matrix, where P stands for the 
%   number of interferograms. It stores the information about how
%   SAR images are connnected to form the interferograms. In each
%   row, it stores the two indexs of SAR images which forms an
%   interferogram, which represents the master image and the slave
%   image respectively. The program will detect the phase closures
%   from this parameter automatically. And if a input interferogram
%   doesn't belong to any phase closure, it will give a warning and
%   that interferogram will be unwrpped separately by a classic 2-D 
%   MCF algorithm. 
INTERFEROGRAMS='ifgs';

%   The coordinate of pixels, which are stored in INPUTFILE. This 
%   parameter is only needed for discrete data to generate TIN. It
%   can be the pixel or the geodetic coordinate. The size of this 
%   parameter should be S*2, where S stands for the number of 
%   pixels in individual interferogram. And in each row, it stores
%   the x and y coordinate of a pixel.
COORDINATE='coor';

%   The output file name. Currently, the unwrapped phases are the only
%   output of this program.
OUTPUTFILE='output_3DMCF_data.mat';

%   Unwrapped phases, which are stored in OUTPUTFILE. It will be the
%   size of WRAPPEDPHASES.
UNWRAPPEDPHASES='uPh';


%================================================================
%   Tile Control. Currently this function can only be used for
%   continuous data, and we simply merge the result of each tile.
%   Therefore, please use this function with caution. And if you
%   don't want to set tile control, setting the NTILER and NTILEC
%   equal one.
%================================================================

%   Tile number in row.
NTILER=1;

%   Tile number in colunm.
NTILEC=1;

%   Overlap between each tile.
OVERLAP=1;




