import numpy as np
import os
import gc
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
import itk
from rollingBall import normalize_local_contrast,rolling_ball_float_background
from skimage import io
import json
import sys

# outputAfterOnlyPoints = False


def replace_nan_with_min(image):
    # Replace NaN values with the minimum value
    image = np.where(np.isfinite(image), image, 0)
    return image

def normalize_local_intensity(imgIn):
    # Setting parameters
    blockRadiusX = 400
    blockRadiusY = 400
    meanFactor = 3.0
    center = True
    stretch = False

    imgIn = normalize_local_contrast(imgIn, blockRadiusX, blockRadiusY, meanFactor, center, stretch)
    # fixedIm = normalize_local_contrast(fixedIm, blockRadiusX, blockRadiusY, meanFactor, center, stretch)

    imgIn = np.array(imgIn)
    imgIn = rolling_ball_float_background(imgIn,50)
    # fixedIm = rolling_ball_float_background(fixedIm,25)


    imgIn = np.array(imgIn)
    blockRadiusX = 50
    blockRadiusY = 50
    meanFactor = 3.0
    center = False
    stretch = True
    imgIn = normalize_local_contrast(imgIn, blockRadiusX, blockRadiusY, meanFactor, center, stretch)
    # fixedIm = normalize_local_contrast(fixedIm, blockRadiusX, blockRadiusY, meanFactor, center, stretch)
    return imgIn

def getCompleteParameterMap(initial_transform_file,initial_affine,numIterations,num_resolutions_penalties,num_resolutions_no_penalties,head,finalGridSpacingInVoxels):
    parameterObjInit = itk.ParameterObject.New() 
    parameterObjInit.ReadParameterFile(initial_transform_file)
    parameterMap = parameterObjInit.GetParameterMap(0)
    parameterMap["TransformParameters"] = tuple([str(num) for num in initial_affine])
    print("affine: " + str(tuple([str(num) for num in initial_affine])))
    nDims = 2 if len(initial_affine)==6 else 3 if len(initial_affine)==12 else None

    parameterMap['FixedImageDimension']=[str(nDims) + "."]
    parameterMap['MovingImageDimension']=[str(nDims) + "."]
    parameterMap['NumberOfParameters']=[str(len(initial_affine)) + "."]
    parameterMap['Origin']=["0."]*nDims
    parameterMap["MaximumNumberOfSamplingAttempts"] = ["500"]
    parameterMap["NumberOfSpatialSamples"] = ["10000"]
    
    parameterObjInit.SetParameterMap(parameterMap)
    parameterObjInit.WriteParameterFile(parameterObjInit,initial_transform_file)

    parameterObj = itk.ParameterObject.New()

    bs_parametermap_penalties = parameterObj.GetDefaultParameterMap('bspline')
    bs_parametermap_penalties["MaximumNumberOfIterations"] = [str(numIterations)]

    bs_parametermap_penalties['NumberOfResolutions'] = [str(num_resolutions_penalties)]

    bs_parametermap_penalties['GridSpacingSchedule'] = [str(int(2**(x-2+num_resolutions_no_penalties))) for x in np.arange(num_resolutions_penalties,0,-1)]
    bs_parametermap_penalties['FinalGridSpacingInPhysicalUnits'] = []
    bs_parametermap_penalties['FinalGridSpacingInVoxels'] = finalGridSpacingInVoxels
    bs_parametermap_penalties['NumberOfSamplesForExactGradient'] = ['100000']
    bs_parametermap_penalties['ImagePyramidSchedule'] = ["1"] *num_resolutions_penalties*nDims
    bs_parametermap_penalties['FixedImagePyramidSchedule'] = bs_parametermap_penalties['ImagePyramidSchedule']
    bs_parametermap_penalties['MovingImagePyramidSchedule'] = bs_parametermap_penalties['ImagePyramidSchedule']
    bs_parametermap_penalties['NumberOfSpatialSamples'] = ['10000']
    bs_parametermap_penalties['NumberOfHistogramBins'] = ['64.0']

    fixed_file_path = os.path.join(head, "points", 'ptsFixed.txt')
    moved_file_path = os.path.join(head, "points", 'ptsMov.txt')
    pointsFilesExist = os.path.isfile(fixed_file_path) and os.path.isfile(moved_file_path)

    if pointsFilesExist:
        print("found pts files")
        original_metric = bs_parametermap_penalties['Metric']
        bs_parametermap_penalties['Registration'] = ['MultiMetricMultiResolutionRegistration']
        bs_parametermap_penalties['Metric'] = [original_metric[0],'CorrespondingPointsEuclideanDistanceMetric','TransformBendingEnergyPenalty']
        bs_parametermap_penalties['Metric0Weight'] = [str(0.1)]
        bs_parametermap_penalties['Metric1Weight'] = [str(0.9)]
        bs_parametermap_penalties['Metric2Weight'] = [str(0.01)]
    else:
        print(fixed_file_path)
        print(moved_file_path)
        raise Exception("points files not found")

    if num_resolutions_penalties>=1:
        parameterObj.AddParameterMap(bs_parametermap_penalties)
    else:
        print("THERE AREN'T PENALTIES RESOLUTIONS")
        print(num_resolutions_penalties)
        
    bs_parametermap_no_penalties = parameterObj.GetDefaultParameterMap('bspline')
    for key in bs_parametermap_penalties.keys():
        bs_parametermap_no_penalties[key] = bs_parametermap_penalties[key]
    bs_parametermap_no_penalties['GridSpacingSchedule'] = [str(int(2**(x-1))) for x in np.arange(num_resolutions_no_penalties,0,-1)]
    bs_parametermap_no_penalties['NumberOfResolutions'] = [str(num_resolutions_no_penalties)]
    bs_parametermap_no_penalties['Metric0Weight'] = [str(1)]
    bs_parametermap_no_penalties['Metric1Weight'] = [str(0)]
    bs_parametermap_no_penalties['Metric2Weight'] = [str(0)]
    bs_parametermap_no_penalties['ImagePyramidSchedule'] = ["1"] * num_resolutions_no_penalties*nDims

    if num_resolutions_no_penalties>=1:
        parameterObj.AddParameterMap(bs_parametermap_no_penalties)
    else:
        print("THERE AREN'T NO PENALTIES RESOLUTIONS")
        print(num_resolutions_no_penalties)
    print(parameterObj)
    return parameterObj,pointsFilesExist,moved_file_path,fixed_file_path,parameterObjInit

def getProcessedImages(imFile,targetChannel,numChannels,adjustContrast):
    im = io.imread(imFile)

    if numChannels > 1:
        channel_dim = im.shape.index(numChannels) if numChannels in im.shape else None
        if channel_dim is not None:
            slices = [slice(im.shape[i]) for i in range(len(im.shape))]
            slices[channel_dim] = int(targetChannel) - 1
            im = im[tuple(slices)]

    im = np.array(im)

    if adjustContrast:
        im = normalize_local_intensity(im)

    im = replace_nan_with_min(im)
    im = itk.GetImageFromArray(im)

    im.SetSpacing([1.0] * im.GetImageDimension())
    return im

def doRegistration(initial_transform_file,fixedIm,movingIm,parameterObj,pointsFilesExist,moved_file_path,fixed_file_path):
    print('running registration')
    
    elastix_object = itk.ElastixRegistrationMethod.New(fixedIm,movingIm)
    elastix_object.SetOutputDirectory('.')
    elastix_object.SetParameterObject(parameterObj)

    if pointsFilesExist:

        elastix_object.SetFixedPointSetFileName(fixed_file_path)
        elastix_object.SetMovingPointSetFileName(moved_file_path)

    for index in range(parameterObj.GetNumberOfParameterMaps()):
        parameter_map = parameterObj.GetParameterMap(index)
        parameterObj.WriteParameterFile(parameter_map, os.path.join(".", "TransformParameters.{0}.txt".format(index)))

    elastix_object.SetInitialTransformParameterFileName(initial_transform_file)
    elastix_object.LogToConsoleOn()
    elastix_object.UpdateLargestPossibleRegion()
    result_transform_parameters = elastix_object.GetTransformParameterObject()
    result_image = elastix_object.GetOutput()
    return result_transform_parameters,result_image

def postProcessRegChannels(fixedIm,result_image):
    
    result_image = result_image-result_image[:].min()
    fixedIm = fixedIm - fixedIm[:].min()

    result_image = (result_image/result_image.max())*(2**15)
    fixedIm = (fixedIm/fixedIm.max())*(2**15)

    result_image = itk.GetImageFromArray(result_image).astype(np.uint16)
    fixedIm = itk.GetImageFromArray(fixedIm).astype(np.uint16)
    composeFilter = itk.ComposeImageFilter.New(fixedIm, result_image)

    result_composed = composeFilter.GetOutput()
    result_composed.Update()
    return result_composed

def updateAndSaveIntermediate(parameterObjInit,result_transform_parameters,head,doInterpolate,nDims,movingIm):
    parameterObjFinal = itk.ParameterObject.New()

    initParameterMap = parameterObjInit.GetParameterMap(0)
    if not doInterpolate:
        initParameterMap["FinalBSplineInterpolationOrder"] = [str(0)]
    parameterObjFinal.AddParameterMap(initParameterMap)

    for index in range(result_transform_parameters.GetNumberOfParameterMaps()):
        parameter_map = result_transform_parameters.GetParameterMap(index)
        parameter_map["InitialTransformParametersFileName"] = ["NoInitialTransform"]
        if movingIm is not None and index == result_transform_parameters.GetNumberOfParameterMaps()-1 and (nDims==3):
            transformix_object = itk.TransformixFilter.New(movingIm)
            transformix_object.SetTransformParameterObject(parameterObjFinal)
            transformix_object.UpdateLargestPossibleRegion()
            result_image = transformix_object.GetOutput()
            result_image = result_image-result_image[:].min()
            result_image = result_image/result_image[:].max()*movingIm[:].max()
            itk.imwrite(itk.GetImageFromArray(result_image.astype(np.uint16)), os.path.join(head,'movingRegChannelOnlyPoints.tif'))
    
        if not doInterpolate:
            parameter_map["FinalBSplineInterpolationOrder"] = [str(0)]
        parameterObjFinal.AddParameterMap(parameter_map)		
    
    for index in range(parameterObjFinal.GetNumberOfParameterMaps()):
        parameter_map = parameterObjFinal.GetParameterMap(index)
        parameterObjFinal.WriteParameterFile(parameter_map, os.path.join(head, "TransformParameters.{0}.txt".format(index)))


    print(parameterObjFinal)
    return parameterObjFinal

def transform_image(image, parameterObjFinal, is_single_channel=True):
    
    image = itk.GetImageFromArray(image)
    transformix_object = itk.TransformixFilter.New(image)
    transformix_object.SetTransformParameterObject(parameterObjFinal)
    transformix_object.UpdateLargestPossibleRegion()
    return transformix_object.GetOutput()

def save_image(image, path):
    itk.imwrite(itk.GetImageFromArray(image.astype(np.uint16)), path)

def process_single_channel_image(image, head, filename, parameterObjFinal):
    transformed_image = transform_image(image, parameterObjFinal)
    save_image(transformed_image, os.path.join(head, "t_" + filename))

def process_multi_channel_image(image, channels, head, parameterObjFinal):
    for k in range(channels):
        print(f'Transforming channel {k+1} of {channels}')
        channel_dim = image.shape.index(channels) if channels in image.shape else None
        print(channel_dim)
        if channel_dim is not None:
            slices = [slice(None)] * len(image.shape)
            slices[channel_dim] = k
            channel_image = image[tuple(slices)]
            transformed_image = transform_image(channel_image, parameterObjFinal, is_single_channel=False)
            save_image(transformed_image, os.path.join(head, f'channel_{k+1}.tif'))

def transformAndSaveOthers(movingImFilename,parameterObjFinal, numChannels, outFile, nDims, fixedImShape):
    head, _ = os.path.split(outFile)
    channels = int(numChannels)

    filename = movingImFilename

    print(filename)
    image = itk.GetImageFromArray(io.imread(filename).astype(np.float32))

    if channels > 1:
        if nDims != 3:
            raise Exception("Channels only supported for 3D images")
        process_multi_channel_image(image, channels, head, parameterObjFinal)
    else:
        if nDims == 3 or len(image.shape) == 2:
            process_single_channel_image(image, head, filename, parameterObjFinal)
        else:
            # Process for each z-slice and save
            transform_result = np.zeros((image.shape[0],) + fixedImShape)
            for k in range(transform_result.shape[0]):
                slice_image = itk.GetImageFromArray(image[k, :, :])
                transformed_slice = transform_image(slice_image, parameterObjFinal)
                transform_result[k, :, :] = transformed_slice[:]
            save_image(transform_result, os.path.join(head, 't_' + filename))


def doIntensityBspline(argsDict):
    print('===========  doIntensityBspline  =================')
    initialAffine = [float(x) for x in argsDict["affineTransform"].split(',')];
    nDims = 2 if len(initialAffine)==6 else 3 if len(initialAffine)==12 else None
    head,tail = os.path.split(argsDict["outputFullFilePath"])


    print('loading images')
    movingIm = getProcessedImages(argsDict["movingImFilename"],argsDict["registrationChannelMoving"],argsDict["movingChannels"],argsDict["preprocessImages"])
    fixedIm = getProcessedImages(argsDict["fixedImFilename"],argsDict["registrationChannelMoving"],argsDict["fixedChannels"],argsDict["preprocessImages"])

    finalGridSpacingInVoxels = [str(float(data["finalGridSpacingxy"]))]*2
    if nDims==3:
        finalGridSpacingInVoxels.append(str(float(data["finalGridSpacingz"])))

    parameterObj,pointsFilesExist,moved_file_path,fixed_file_path,parameterObjInit = getCompleteParameterMap(argsDict["initialTransformFilePath"],initialAffine,argsDict["iterations"],argsDict["resolutionsPenalties"],argsDict["resolutionsNoPenalties"],head,finalGridSpacingInVoxels)

    result_transform_parameters,result_image = doRegistration(argsDict["initialTransformFilePath"],fixedIm,movingIm,parameterObj,pointsFilesExist,moved_file_path,fixed_file_path)
    gc.collect()
    
    result_composed = postProcessRegChannels(fixedIm,result_image)
    itk.imwrite(result_composed, argsDict["outputFullFilePath"])
    fixedImShape = fixedIm.shape # Save for later
    del result_composed, fixedIm, result_image  # Free memory
    gc.collect()

    # if not outputAfterOnlyPoints:
    #     movingIm = None
    
    parameterObjFinal = updateAndSaveIntermediate(parameterObjInit,result_transform_parameters,head,argsDict["interpolateTransform"],nDims,movingIm)
    
    transformAndSaveOthers(argsDict["movingImFilename"],parameterObjFinal,argsDict["toTransformChannels"],argsDict["outputFullFilePath"],nDims,fixedImShape)

    return


if __name__ == "__main__":
    # jsonFile = sys.argv[1]
    jsonFile = r'C:\Users\jpv88\Documents\GitHub\HCRprocess\internal\data.json'
    with open(jsonFile) as f:
        data = json.load(f)
        print(data)
        doIntensityBspline(data)
        print('done')