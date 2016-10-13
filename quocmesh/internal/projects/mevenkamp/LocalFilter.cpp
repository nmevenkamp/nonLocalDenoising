#include <localFilter.h>
#include <statistics.h>
#include <multiStreambuf.h>

typedef double RType;
typedef aol::FullMatrix<RType> MType;
typedef qc::ScalarArray<RType, qc::QC_2D> PicType;
typedef qc::MultiArray<RType, qc::QC_2D, 3> ColoredPicType;

typedef RType RealType;
typedef MType MatrixType;
typedef PicType PictureType;
typedef ColoredPicType ColoredPictureType;


int main ( int argc, char** argv ) {
  try {
    aol::ParameterParser parser ( argc, argv, "localFilter.par" );
    const bool quietMode = parser.checkAndGetBool ( "quietMode" );
    
    const std::string srcPath = parser.getString ( "sourcePath" );
    const std::string outputDir = parser.getString ( "outputDir" );
    const std::string groundTruthPath = parser.getStringOrDefault ( "groundTruthPath", "" );
    bool groundTruthAvailable = groundTruthPath != "" && aol::fileExists ( groundTruthPath );
    
    const int filterType = parser.getIntOrDefault ( "filterType", 0 );
    
    const int noiseType = im::resolveIdentifier<im::NOISE_TYPE> ( parser.getStringOrDefault ( "noiseType", im::NOISE_TYPE::getIdentifier ( im::NOISE_TYPE::GAUSSIAN ) ) );
    
    const short blockSize = parser.getIntOrDefault ( "blockSize", 5 );
    if ( blockSize <= 0 )
      throw aol::Exception ( "Block size must be greater than or equal to 1!", __FILE__, __LINE__ );
    
    
    PictureType noisy ( srcPath.c_str ( ) ), groundTruth, estimate ( noisy, aol::STRUCT_COPY );
    if ( groundTruthAvailable ) groundTruth.load ( groundTruthPath.c_str ( ) );
    groundTruthAvailable = groundTruth.getNumX ( ) == noisy.getNumX ( ) && groundTruth.getNumY ( ) == noisy.getNumY ( );
    
    aol::ProgressBar<> progressBar ( "Denoising", std::cerr );
    im::LocalFilter<RealType> localFilter;
    im::LocalFilterOptions<RealType> options ( "", quietMode );
    if ( groundTruthAvailable ) options.groundTruth = &groundTruth;
    options.progressBar = &progressBar;
    options.filterType = filterType;
    options.noiseType = noiseType;
    options.alpha = parser.getDoubleOrDefault ( "alpha", 0.0 );
    options.sigma = parser.getDoubleOrDefault ( "sigma", 0.0 );
    options.mu = parser.getDoubleOrDefault ( "mu", 0.0 );
    options.blockSize.setAll ( blockSize );
    
    std::string dirPath = aol::strprintf ( "%s/%s_%s", outputDir.c_str ( ), aol::getBaseFileName ( srcPath ).c_str ( ), localFilter.getMethodIdentifier ( options ).c_str ( ) );
    aol::makeDirectory ( dirPath.c_str ( ) );
    options.outputDir = dirPath;
    aol::AdditionalOutputToFile addOut ( aol::strprintf ( "%s/log.txt", dirPath.c_str ( ) ).c_str ( ) );
    
    
    aol::StopWatch watch;
    watch.start ( );
    localFilter.apply ( options, noisy, estimate );
    watch.stop ( );
    if ( !quietMode ) watch.printReport ( std::cerr );
    
    noisy.save ( aol::strprintf ( "%s/noisy%s", dirPath.c_str ( ), getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    estimate.save ( aol::strprintf ( "%s/estimate%s", dirPath.c_str ( ), getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    if ( groundTruthAvailable )
      groundTruth.save ( aol::strprintf ( "%s/groundTruth%s", dirPath.c_str ( ), getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    
    
    if ( filterType == 0 || parser.checkAndGetBool ( "cropBordersBeforeAnalysis" ) ) {
      aol::Vec<2, int> cropStart, cropSize;
      cropStart[0] = ( blockSize - 1 ) / 2;
      cropStart[1] = ( blockSize - 1 ) / 2;
      cropSize[0] = noisy.getNumX ( ) - ( blockSize - 1 );
      cropSize[1] = noisy.getNumY ( ) - ( blockSize - 1 );
      noisy.crop ( cropStart, cropSize );
      groundTruth.crop ( cropStart, cropSize );
      estimate.crop ( cropStart, cropSize );
    }
    
    
    if ( !quietMode && groundTruthAvailable ) {
      std::cerr << "MSE:\t" << aol::MSE<RealType> ( groundTruth, estimate ) << std::endl;
      if ( noiseType == im::NOISE_TYPE::GAUSSIAN ) std::cerr << "PSNR:\t" << aol::PSNR<RealType> ( groundTruth, estimate, 255 ) << std::endl;
      else if ( noiseType == im::NOISE_TYPE::POISSON || noiseType == im::NOISE_TYPE::CMOS ) std::cerr << "PSNR:\t" << aol::PSNR<RealType> ( groundTruth, estimate ) << std::endl;
    }
    
    if ( noiseType == im::NOISE_TYPE::GAUSSIAN ) {
      saveMethodNoise ( noisy, estimate, aol::strprintf ( "%s/methodNoise%s", dirPath.c_str ( ), getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ) );
      saveInverseImageGaussian<RealType, PictureType> ( noisy, estimate, localFilter.sigma ( ), dirPath );
    } else if ( noiseType == im::NOISE_TYPE::POISSON ) {
      saveInverseImagePoisson<RealType, PictureType> ( noisy, estimate, aol::strprintf ( "%s/inverseImage%s", dirPath.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ) );
      saveInverseHistogramPoisson<RealType, PictureType> ( noisy, estimate, aol::strprintf ( "%s/inverseHistogram", dirPath.c_str ( ) ).c_str ( ) );
      if ( !quietMode && groundTruthAvailable ) {
        saveInverseImagePoisson<RealType, PictureType> ( noisy, groundTruth, aol::strprintf ( "%s/inverseImage_gt%s", dirPath.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ) );
        saveInverseHistogramPoisson<RealType, PictureType> ( noisy, groundTruth, aol::strprintf ( "%s/inverseHistogram_gt", dirPath.c_str ( ) ).c_str ( ) );
      }
    } else if ( noiseType == im::NOISE_TYPE::CMOS ) {
      saveInverseImageMixedPoissonGaussian<RealType, PictureType> ( noisy, estimate, localFilter.alpha ( ), localFilter.mu ( ), localFilter.sigma ( ), aol::strprintf ( "%s/inverseImage%s", dirPath.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ) );
      if ( !quietMode && groundTruthAvailable )
        saveInverseImageMixedPoissonGaussian<RealType, PictureType> ( noisy, groundTruth, localFilter.alpha ( ), localFilter.mu ( ), localFilter.sigma ( ), aol::strprintf ( "%s/inverseImage_gt%s", dirPath.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ) );
    }
  } catch ( aol::Exception& ex ) {
    ex.dump();
  }
}