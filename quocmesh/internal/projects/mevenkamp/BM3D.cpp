#include <bm3dFilter.h>
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
    aol::ParameterParser parser ( argc, argv, "bm3d.par" );
    const bool quietMode = parser.checkAndGetBool ( "quietMode" );
    
    const std::string srcPath = parser.getStringExpandTilde ( "sourcePath" );
    const std::string outputDir = parser.getStringExpandTilde ( "outputDir" );
    const std::string groundTruthPath = parser.getStringOrDefault ( "groundTruthPath", "" );
    bool groundTruthAvailable = groundTruthPath != "" && aol::fileExists ( groundTruthPath );

    const int noiseType = im::resolveIdentifier<im::NOISE_TYPE> ( parser.getStringOrDefault ( "noiseType", im::NOISE_TYPE::getIdentifier ( im::NOISE_TYPE::GAUSSIAN ) ) );
    
    const int poissonNoiseAdaptation = im::resolveIdentifier<im::POISSON_NOISE_ADAPTATION> ( parser.getStringOrDefault ( "poissonNoiseAdaptation",
      im::POISSON_NOISE_ADAPTATION::getIdentifier ( im::POISSON_NOISE_ADAPTATION::ANSCOMBE ) ) );
    
    const int profile = im::resolveIdentifier<im::BM3D_PROFILE> ( parser.getStringOrDefault ( "profile", im::BM3D_PROFILE::getIdentifier ( im::BM3D_PROFILE::NP ) ) );
    
    short blockSize = parser.getIntOrDefault ( "blockSize", 8 );
 
    const int similaritySearchMethod = im::resolveIdentifier<im::BM3D_SIMILARITYSEARCH_METHOD> ( parser.getStringOrDefault ( "similaritySearchMethod",
      im::BM3D_SIMILARITYSEARCH_METHOD::getIdentifier ( im::BM3D_SIMILARITYSEARCH_METHOD::LOCAL ) ) );
    
    PictureType noisy ( srcPath.c_str ( ) ), groundTruth, estimate;
    if ( groundTruthAvailable ) groundTruth.load ( groundTruthPath.c_str ( ) );
    groundTruthAvailable = groundTruth.getNumX ( ) == noisy.getNumX ( ) && groundTruth.getNumY ( ) == noisy.getNumY ( );
    
    aol::ProgressBar<> progressBar ( "Denoising", std::cerr );
    im::BM3DFilter<RealType> bm3dFilter;
    im::BM3DOptions<RealType> options ( "", quietMode );
    if ( groundTruthAvailable ) options.groundTruth = &groundTruth;
    options.progressBar = &progressBar;
    options.noiseType = noiseType;
    options.alpha = parser.getDoubleOrDefault ( "alpha", 0.0 );
    options.sigma = parser.getDoubleOrDefault ( "sigma", 0.0 );
    options.mu = parser.getDoubleOrDefault ( "mu", 0.0 );
    if ( options.noiseType == im::NOISE_TYPE::GAUSSIAN || options.noiseType == im::NOISE_TYPE::CMOS ) {
      options.levelSetNumSamples = parser.getIntOrDefault ( "levelSetNumSamples", 100 );
      options.levelSetNumPointsPerSample = parser.getIntOrDefault ( "levelSetNumPointsPerSample", 64 );
      options.levelSetOverSamplingFactor = parser.getIntOrDefault ( "levelSetOverSamplingFactor", 10 );
    }
    if ( noiseType == im::NOISE_TYPE::POISSON ) options.poissonNoiseAdaptation = poissonNoiseAdaptation;
    options.profile = profile;
    options.blockSize.setAll ( blockSize );
    options.similaritySearchMethod = similaritySearchMethod;
    options.denoiseInTransformDomain = parser.getBoolOrDefault ( "denoiseInTransformDomain", true );
    
    std::cerr << bm3dFilter.getMethodIdentifier ( options ) << std::endl;
    std::string dirPath = aol::strprintf ( "%s/%s_%s", outputDir.c_str ( ), aol::getBaseFileName ( srcPath ).c_str ( ), bm3dFilter.getMethodIdentifier ( options ).c_str ( ) );
    aol::makeDirectory ( dirPath.c_str ( ) );
    options.outputDir = dirPath;
    aol::AdditionalOutputToFile addOut ( aol::strprintf ( "%s/log.txt", dirPath.c_str ( ) ).c_str ( ) );

    bm3dFilter.setCatchCtrlC ( true );
    
    if ( !quietMode ) std::cerr << "Denoising image.." << std::endl;
    aol::StopWatch watch;
    watch.start ( );

    bm3dFilter.apply ( options, noisy, estimate );
    
    watch.stop ( );
    if ( !quietMode ) watch.printReport ( std::cerr );
    
    noisy.save ( aol::strprintf ( "%s/noisy%s", dirPath.c_str ( ), getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    estimate.save ( aol::strprintf ( "%s/estimate%s", dirPath.c_str ( ), getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    if ( groundTruthAvailable )
      groundTruth.save ( aol::strprintf ( "%s/groundTruth%s", dirPath.c_str ( ), getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ), qc::PGM_DOUBLE_BINARY );

    // Crop images (border is not denoised by NLM currently)
    if ( parser.checkAndGetBool ( "cropBordersBeforeAnalysis" ) ) {
      blockSize = aol::Max<int> ( 1, bm3dFilter.blockSize ( )[0] );
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
      saveInverseImageGaussian<RealType, PictureType> ( noisy, estimate, bm3dFilter.sigma ( ), dirPath );
    } else if ( noiseType == im::NOISE_TYPE::POISSON ) {
      saveInverseImagePoisson<RealType, PictureType> ( noisy, estimate, aol::strprintf ( "%s/inverseImage%s", dirPath.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ) );
      saveInverseHistogramPoisson<RealType, PictureType> ( noisy, estimate, aol::strprintf ( "%s/inverseHistogram", dirPath.c_str ( ) ).c_str ( ) );
      if ( !quietMode && groundTruthAvailable ) {
        saveInverseImagePoisson<RealType, PictureType> ( noisy, groundTruth, aol::strprintf ( "%s/inverseImage_gt%s", dirPath.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ) );
        saveInverseHistogramPoisson<RealType, PictureType> ( noisy, groundTruth, aol::strprintf ( "%s/inverseHistogram_gt", dirPath.c_str ( ) ).c_str ( ) );
      }
    } else if ( noiseType == im::NOISE_TYPE::CMOS ) {
      saveInverseImageMixedPoissonGaussian<RealType, PictureType> ( noisy, estimate, bm3dFilter.alpha ( ), bm3dFilter.mu ( ), bm3dFilter.sigma ( ), aol::strprintf ( "%s/inverseImage%s", dirPath.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ) );
      if ( !quietMode && groundTruthAvailable )
        saveInverseImageMixedPoissonGaussian<RealType, PictureType> ( noisy, groundTruth, bm3dFilter.alpha ( ), bm3dFilter.mu ( ), bm3dFilter.sigma ( ), aol::strprintf ( "%s/inverseImage_gt%s", dirPath.c_str ( ), qc::getDefaultArraySuffix ( qc::QC_2D ) ).c_str ( ) );
    }
  } catch ( aol::Exception& ex ) {
    ex.dump();
  }
}