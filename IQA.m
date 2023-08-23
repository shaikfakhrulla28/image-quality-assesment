function varargout = IQA(varargin)

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IQA_OpeningFcn, ...
                   'gui_OutputFcn',  @IQA_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before IQA is made visible.
function IQA_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

ResetButton_Callback(hObject, eventdata, handles)

% --- Outputs from this function are returned to the command line.
function varargout = IQA_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in BrowseImage.
function BrowseImage_Callback(hObject, eventdata, handles)

ResetButton_Callback(hObject, eventdata, handles);

global image;
[filename pathname] = uigetfile({'*.jpg';'*.bmp';'*.tif';'*.png'},'File Selector');
x = strcat(pathname, filename);
image=imread(x);
axes(handles.axes1);
imshow(image);

% --- Executes on button press in AddNoise.
function AddNoise_Callback(hObject, eventdata, handles)

global image;
global addnoisyimage;
global mean;
global variance;
global AdditiveNoiseMenu;
if (strcmp(AdditiveNoiseMenu, 'Gaussian'))
    addnoisyimage = imnoise(image, 'Gaussian', mean, variance);
elseif (strcmp(AdditiveNoiseMenu, 'Poisson'))
    addnoisyimage = imnoise(image, 'Poisson');
    elseif (strcmp(AdditiveNoiseMenu, 'Select Additive Noise Type'))
    addnoisyimage = image;
end
axes(handles.axes2);
imshow(addnoisyimage);

% --- Executes on button press in MultiNoise.
function MultiNoise_Callback(hObject, eventdata, handles)
global noisedensity;
global variance_multi;
global image;
global multinoisyimage;
global MultiplicativeNoiseMenu;
if (strcmp(MultiplicativeNoiseMenu, 'Salt & Pepper'))
    multinoisyimage = imnoise(image, 'salt & pepper', noisedensity);
elseif (strcmp(MultiplicativeNoiseMenu, 'Speckle'))
    multinoisyimage = imnoise(image, 'speckle', variance_multi);
elseif (strcmp(MultiplicativeNoiseMenu, 'Select Multiplicative Noise'))
    multinoisyimage = image;
end
axes(handles.axes3);
imshow(multinoisyimage);

% --- Executes on button press in CheckPSNR.
function CheckPSNR_Callback(hObject, eventdata, handles)

global addnoisyimage;
global multinoisyimage;
global image;
global s;
global u;
global justforcontrol;

if (get(hObject, 'Value') == get(hObject,'Max'))
    justforcontrol=1;
    s=psnr(addnoisyimage, image);
    u=psnr(multinoisyimage, image);
else
    justforcontrol=0;
    s='--';
    u='--';
end

% Hint: get(hObject,'Value') returns toggle state of CheckPSNR
function s = psnr(addnoisyimage, image)

if(ndims(addnoisyimage)==3)
    addnoisyimage = rgb2gray(addnoisyimage);
end

if(ndims(image)==3)
    image = rgb2gray(image);
end

addnoisyimage=double(addnoisyimage);
image=double(image);

[m,n] = size(addnoisyimage);

peak=255*255*m*n;

noise  = addnoisyimage - image;
nostotal = sum(sum(noise.*noise));

if nostotal == 0
    s = 'INF'; %% INF. clean image
else
    s = 10 * log10(peak./nostotal);
end

% --- Executes on button press in CheckSSIM.
function CheckSSIM_Callback(hObject, eventdata, handles)

global addnoisyimage;
global multinoisyimage;
global image;
global t;
global v;
global justforcontrol2;
K = [0.05 0.05];
window = ones(8);
L = 100;
Z = [0.01 0.03];
if (get(hObject, 'Value') == get(hObject,'Max'))
    justforcontrol2=1;
    t=ssim(addnoisyimage, image, Z, window, L);
    v=ssim(multinoisyimage, image, Z, window, L);
else
    justforcontrol2=0;
    t='--';
    v='--';
end

% Hint: get(hObject,'Value') returns toggle state of CheckSSIM
function [mssim] = ssim(img1, img2, Z, window, L)

if(ndims(img1)==3)
    img1=rgb2gray(img1);
end
if(ndims(img2)==3)
    img2=rgb2gray(img2);
end

[rows,cols]=size(img2);
img1=imresize(img1,[rows cols]);

if (nargin < 2 || nargin > 5)
   mssim = -Inf;
   ssim_map = -Inf;
   return;
end

if (size(img1) ~= size(img2))
   mssim = -Inf;
   ssim_map = -Inf;
   return;
end

[M N] = size(img1);

if (nargin == 2)
   if ((M < 11) || (N < 11))
	   mssim = -Inf;
	   ssim_map = -Inf;
      return
   end
   window = fspecial('gaussian', 11, 1.5);	%
   Z(1) = 0.01;					% default settings
   Z(2) = 0.03;			
   L = 255;                          
end

if (nargin == 3)
   if ((M < 11) || (N < 11))
	   mssim = -Inf;
	   ssim_map = -Inf;
      return
   end
   window = fspecial('gaussian', 11, 1.5);
   L = 255;
   if (length(Z) == 2)
      if (Z(1) < 0 || Z(2) < 0)
		   mssim = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   mssim = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

if (nargin == 4)
   [H W] = size(window);
   if ((H*W) < 4 || (H > M) || (W > N))
	   mssim = -Inf;
	   ssim_map = -Inf;
      return
   end
   L = 255;
   if (length(Z) == 2)
      if (Z(1) < 0 || Z(2) < 0)
		   mssim = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   mssim = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

if (nargin == 5)
   [H W] = size(window);
   if ((H*W) < 4 || (H > M) || (W > N))
	   mssim = -Inf;
	   ssim_map = -Inf;
      return
   end
   if (length(Z) == 2)
      if (Z(1) < 0 || Z(2) < 0)
		   mssim = -Inf;
   		ssim_map = -Inf;
	   	return;
      end
   else
	   mssim = -Inf;
   	ssim_map = -Inf;
	   return;
   end
end

img1 = double(img1);
img2 = double(img2);

% automatic downsampling
f = max(1,round(min(M,N)/256));
%downsampling by f
%use a simple low-pass filter 
if(f>1)
    lpf = ones(f,f);
    lpf = lpf/sum(lpf(:));
    img1 = imfilter(img1,lpf,'symmetric','same');
    img2 = imfilter(img2,lpf,'symmetric','same');

    img1 = img1(1:f:end,1:f:end);
    img2 = img2(1:f:end,1:f:end);
end

C1 = (Z(1)*L)^2;
C2 = (Z(2)*L)^2;
window = window/sum(sum(window));

mu1   = filter2(window, img1, 'valid');
mu2   = filter2(window, img2, 'valid');
mu1_sq = mu1.*mu1;
mu2_sq = mu2.*mu2;
mu1_mu2 = mu1.*mu2;
sigma1_sq = filter2(window, img1.*img1, 'valid') - mu1_sq;
sigma2_sq = filter2(window, img2.*img2, 'valid') - mu2_sq;
sigma12 = filter2(window, img1.*img2, 'valid') - mu1_mu2;

if (C1 > 0 && C2 > 0)
   ssim_map = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2));
else
   numerator1 = 2*mu1_mu2 + C1;
   numerator2 = 2*sigma12 + C2;
	denominator1 = mu1_sq + mu2_sq + C1;
   denominator2 = sigma1_sq + sigma2_sq + C2;
   ssim_map = ones(size(mu1));
   index = (denominator1.*denominator2 > 0);
   ssim_map(index) = (numerator1(index).*numerator2(index))./(denominator1(index).*denominator2(index));
   index = (denominator1 ~= 0) & (denominator2 == 0);
   ssim_map(index) = numerator1(index)./denominator1(index);
end

mssim = mean2(ssim_map);

function edit1_Callback(hObject, eventdata, handles)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in AddNoiseResults.
function AddNoiseResults_Callback(hObject, eventdata, handles)

global u;
global v;
global s;
global t;
global justforcontrol;
global justforcontrol2;

if justforcontrol==1;
    CheckPSNR_Callback(hObject, eventdata, handles)
end
if justforcontrol2==1;
    CheckSSIM_Callback(hObject, eventdata, handles)
end
set(handles.PSNR_Add, 'string', s);
set(handles.SSIM_Add, 'string', t);
set(handles.PSNR_Multi, 'string', u);
set(handles.SSIM_Multi, 'string', v);

% --- Executes on button press in ResetButton.
function ResetButton_Callback(hObject, eventdata, handles)

global s;
global t
global u;
global v;
global q;
global justforcontrol;
global justforcontrol2;
global addnoisyimage;
global multinoisyimage;
global image;
global mean;
global variance;
global noisedensity;
global variance_multi;
global AdditiveNoiseMenu;
global MultiplicativeNoiseMenu;
t='--';
s='--';
u='--';
v='--';
q=0;
justforcontrol=0;
justforcontrol2=0;
image = ones(600,400);
addnoisyimage = image;
multinoisyimage = image;
axes(handles.axes1);
imshow(image);
axes(handles.axes2);
imshow(addnoisyimage);
axes(handles.axes3);
imshow(multinoisyimage);
set(handles.CheckPSNR, 'Value', q);
set(handles.CheckSSIM, 'Value', q);
set(handles.PSNR_Add, 'string', s);
set(handles.SSIM_Add, 'string', t);
set(handles.PSNR_Multi, 'string', u);
set(handles.SSIM_Multi, 'string', v);
set(handles.MeanValue, 'string', mean);
set(handles.VarianceValue, 'string', variance);
set(handles.NoiseDensityValue, 'string', noisedensity);
set(handles.VarianceValue_Multi, 'string', variance_multi);
if (strcmp(AdditiveNoiseMenu, 'Poisson'))
     mean = 'N/A';
     variance = 'N/A';
     set(handles.MeanValue, 'string', mean);
     set(handles.VarianceValue, 'string', variance);
elseif (strcmp(AdditiveNoiseMenu, 'Gaussian'))
    mean=0;
    variance=0.01;
    set(handles.MeanValue, 'string', mean);
    set(handles.VarianceValue, 'string', variance);
elseif (strcmp(AdditiveNoiseMenu, 'Select Additive Noise Type'))
    mean=0;
    variance=0;
    set(handles.MeanValue, 'string', mean);
    set(handles.VarianceValue, 'string', variance);
end
if (strcmp(MultiplicativeNoiseMenu, 'Speckle'))
     noisedensity = 'N/A';
     variance_multi = 0.04;
     set(handles.NoiseDensityValue, 'string', noisedensity);
     set(handles.VarianceValue_Multi, 'string', variance_multi);
elseif (strcmp(MultiplicativeNoiseMenu, 'Salt & Pepper'))
     noisedensity = 0.05;
     variance_multi = 'N/A';
     set(handles.NoiseDensityValue, 'string', noisedensity);
     set(handles.VarianceValue_Multi, 'string', variance_multi);
elseif (strcmp(MultiplicativeNoiseMenu, 'Select Multiplicative Noise'))
    noisedensity = 0;
    variance_multi = 0;
    set(handles.NoiseDensityValue, 'string', noisedensity);
    set(handles.VarianceValue_Multi, 'string', variance_multi);
end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over PSNR_Add.
function PSNR_Add_ButtonDownFcn(hObject, eventdata, handles)

function MeanValue_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of MeanValue as text
%        str2double(get(hObject,'String')) returns contents of MeanValue as a double
global mean;
mean = str2double(get(handles.MeanValue,'string'));

% --- Executes during object creation, after setting all properties.
function MeanValue_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function VarianceValue_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of VarianceValue as text
%        str2double(get(hObject,'String')) returns contents of VarianceValue as a double
global variance;
variance = str2double(get(handles.VarianceValue,'string'));

% --- Executes during object creation, after setting all properties.
function VarianceValue_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function NoiseDensityValue_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of NoiseDensityValue as text
%        str2double(get(hObject,'String')) returns contents of NoiseDensityValue as a double
global noisedensity;
noisedensity = str2double(get(handles.NoiseDensityValue,'string'));

% --- Executes during object creation, after setting all properties.
function NoiseDensityValue_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in AdditiveNoiseMenu.
function AdditiveNoiseMenu_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns AdditiveNoiseMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from AdditiveNoiseMenu
global AdditiveNoiseMenu;
global mean;
global variance;
contents = cellstr(get(hObject,'String'));
AdditiveNoiseMenu = contents{get(hObject,'Value')};
if (strcmp(AdditiveNoiseMenu, 'Poisson'))
    mean = 'N/A';
    variance = 'N/A';
    set(handles.MeanValue, 'string', mean);
    set(handles.VarianceValue, 'string', variance);
elseif (strcmp(AdditiveNoiseMenu, 'Gaussian'))
    mean=0;
    variance=0.01;
    set(handles.MeanValue, 'string', mean);
    set(handles.VarianceValue, 'string', variance);
end

% --- Executes during object creation, after setting all properties.
function AdditiveNoiseMenu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in MultiplicativeNoiseMenu.
function MultiplicativeNoiseMenu_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns MultiplicativeNoiseMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MultiplicativeNoiseMenu
global MultiplicativeNoiseMenu;
global noisedensity;
global variance_multi;
contents = cellstr(get(hObject,'String'));
MultiplicativeNoiseMenu = contents{get(hObject,'Value')};
if (strcmp(MultiplicativeNoiseMenu, 'Speckle'))
    noisedensity = 'N/A';
    variance_multi = 0.04;
    set(handles.NoiseDensityValue, 'string', noisedensity);
    set(handles.VarianceValue_Multi, 'string', variance_multi);
elseif (strcmp(MultiplicativeNoiseMenu, 'Salt & Pepper'))
    noisedensity = 0.05;
    variance_multi = 'N/A';
    set(handles.NoiseDensityValue, 'string', noisedensity);
    set(handles.VarianceValue_Multi, 'string', variance_multi);
end
% --- Executes during object creation, after setting all properties.
function MultiplicativeNoiseMenu_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function VarianceValue_Multi_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of VarianceValue_Multi as text
%        str2double(get(hObject,'String')) returns contents of VarianceValue_Multi as a double
global variance_multi;
variance_multi = str2double(get(handles.VarianceValue_Multi,'string'));

% --- Executes during object creation, after setting all properties.
function VarianceValue_Multi_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
