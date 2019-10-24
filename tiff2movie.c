// Reads tiffs and converts them to a movie
// Uses shell ordering
// 20.6.2017
// Jurij Kotar
// Compile with gcc -o tiff2movie tiff2movie.c -O3 -Wall -ltiff -g
// add the frame rate of your videos at the variable 'double time_interval' 

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <tiffio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <stdbool.h>
#include <libgen.h>

//
// Common camera defines
//
#define CAMERA_MOVIE_MAGIC 0x496d6554 // TemI
#define CAMERA_MOVIE_VERSION 1
#define CAMERA_TYPE_IIDC	1
#define CAMERA_TYPE_ANDOR	2
#define CAMERA_TYPE_XIMEA	3
#define CAMERA_PIXELMODE_MONO_8		1
#define CAMERA_PIXELMODE_MONO_16BE	2 // Big endian
#define CAMERA_PIXELMODE_MONO_16LE	3 // Little endian

//
// IIDC defines
//
#define IIDC_MOVIE_HEADER_LENGTH 172
// Feature modes
#define IIDC_FEATURE_MODE_OFF ( 1<<0 )
#define IIDC_FEATURE_MODE_RELATIVE ( 1<<1 )
#define IIDC_FEATURE_MODE_ABSOLUTE ( 1<<2 )
#define IIDC_FEATURE_MODE_AUTO ( 1<<3 )
#define IIDC_FEATURE_MODE_ONEPUSH ( 1<<4 )
#define IIDC_FEATURE_MODE_ADVANCED ( 1<<5 )
// Trigger
#define IIDC_TRIGGER_INTERNAL  -1
#define IIDC_TRIGGER_EXTERNAL0 0
#define IIDC_TRIGGER_EXTERNAL1 1
#define IIDC_TRIGGER_EXTERNAL15 7

//
// Andor defines
//
#define ANDOR_MOVIE_HEADER_LENGTH 128
// VS Speeds
#define ANDOR_VALUE_VS_SPEED_MIN 4
#define ANDOR_VALUE_VS_SPEED_MAX 0
#define ANDOR_VALUE_VS_SPEED_0_3 0
#define ANDOR_VALUE_VS_SPEED_0_5 1
#define ANDOR_VALUE_VS_SPEED_0_9 2
#define ANDOR_VALUE_VS_SPEED_1_7 3
#define ANDOR_VALUE_VS_SPEED_3_3 4
// VS Amplitudes
#define ANDOR_VALUE_VS_AMPLITUDE_MIN 0
#define ANDOR_VALUE_VS_AMPLITUDE_MAX 4
#define ANDOR_VALUE_VS_AMPLITUDE_0 0
#define ANDOR_VALUE_VS_AMPLITUDE_1 1
#define ANDOR_VALUE_VS_AMPLITUDE_2 2
#define ANDOR_VALUE_VS_AMPLITUDE_3 3
#define ANDOR_VALUE_VS_AMPLITUDE_4 4
// Shutter
#define ANDOR_VALUE_SHUTTER_AUTO 0
#define ANDOR_VALUE_SHUTTER_OPEN 1
#define ANDOR_VALUE_SHUTTER_CLOSE 2
// Cooler
#define ANDOR_VALUE_COOLER_OFF 0
#define ANDOR_VALUE_COOLER_ON 1
// Cooler mode
#define ANDOR_VALUE_COOLER_MODE_RETURN 0
#define ANDOR_VALUE_COOLER_MODE_MAINTAIN 1
// Fan
#define ANDOR_VALUE_FAN_FULL 0
#define ANDOR_VALUE_FAN_LOW 1
#define ANDOR_VALUE_FAN_OFF 2
// ADC
#define ANDOR_VALUE_ADC_14BIT 0
#define ANDOR_VALUE_ADC_16BIT 1
// Amplifier
#define ANDOR_VALUE_AMPLIFIER_EM 0
#define ANDOR_VALUE_AMPLIFIER_CON 1
// Preamp gain
#define ANDOR_VALUE_PREAMP_GAIN_1_0 0
#define ANDOR_VALUE_PREAMP_GAIN_2_4 1
#define ANDOR_VALUE_PREAMP_GAIN_5_1 2
// Trigger
#define ANDOR_VALUE_TRIGGER_INTERNAL  0
#define ANDOR_VALUE_TRIGGER_EXTERNAL 1
#define ANDOR_VALUE_TRIGGER_FAST_EXTERNAL -1 // Combination of external and SetFastExtTrigger
#define ANDOR_VALUE_TRIGGER_EXTERNAL_START 6
#define ANDOR_VALUE_TRIGGER_EXTERNAL_EXPOSURE  7
#define ANDOR_VALUE_TRIGGER_SOFTWARE  10

//
// Ximea defines
//
#define XIMEA_MOVIE_HEADER_LENGTH 64

//
// Common camera frame struct
//
struct camera_save_struct {
	//
	// Common stuff
	//
	uint32_t magic; // 'AndO'
	uint32_t version;
	uint32_t type; // Camera type
	uint32_t pixelmode; // Pixel mode
	uint32_t length_header; // Header data in bytes ( Everything except image data )
	uint32_t length_data; // Total data length in bytes;
};

//
// IIDC movie frame struct
//
union iidc_save_feature_value {
	uint32_t value;
	float absvalue;
};

struct iidc_save_struct {
	//
	// Common stuff
	//
	uint32_t magic; // 'TemI'
	uint32_t version;
	uint32_t type; // Camera type
	uint32_t pixelmode; // Pixel mode
	uint32_t length_header; // Header data in bytes ( Everything except image data )
	uint32_t length_data; // Total data length in bytes;

	//
	// Camera specific stuff
	//
	// Camera properties
	uint64_t i_guid;
	uint32_t i_vendor_id;
	uint32_t i_model_id;

	// Frame properties
	uint32_t i_video_mode;
	uint32_t i_color_coding;

	uint64_t i_timestamp; // microseconds

	uint32_t i_size_x_max; // Sensor size
	uint32_t i_size_y_max;
	uint32_t i_size_x; // Selected region
	uint32_t i_size_y;
	uint32_t i_pos_x;
	uint32_t i_pos_y;

	uint32_t i_pixnum; // Number of pixels
	uint32_t i_stride; // Number of bytes per image line
	uint32_t i_data_depth;  // Number of bits per pixel.

	uint32_t i_image_bytes; // Number of bytes used for the image (image data only, no padding)
	uint64_t i_total_bytes; // Total size of the frame buffer in bytes. May include packet multiple padding and intentional padding (vendor specific)

	// Features
	uint32_t i_brightness_mode; // Current mode
	union iidc_save_feature_value i_brightness; // Can be also float if mode is IIDC_FEATURE_MODE_ABSOLUTE (1<<2)

	uint32_t i_exposure_mode;
	union iidc_save_feature_value i_exposure;

	uint32_t i_gamma_mode;
	union iidc_save_feature_value i_gamma;

	uint32_t i_shutter_mode;
	union iidc_save_feature_value i_shutter;

	uint32_t i_gain_mode;
	union iidc_save_feature_value i_gain;

	uint32_t i_temperature_mode;
	union iidc_save_feature_value i_temperature;

	uint32_t i_trigger_delay_mode;
	union iidc_save_feature_value i_trigger_delay;

	int32_t i_trigger_mode;

	// Advanced features
	uint32_t i_avt_channel_balance_mode;
	int32_t i_avt_channel_balance;

	// Image data
	uint8_t *data;
} __attribute__((__packed__));

//
// Andor movie frame struct
//
struct andor_save_struct {
	//
	// Common stuff
	//
	uint32_t magic; // 'TemI'
	uint32_t version;
	uint32_t type; // Camera type
	uint32_t pixelmode; // Pixel mode
	uint32_t length_header; // Header data in bytes ( Everything except image data )
	uint32_t length_data; // Total data length in bytes;

	//
	// Camera specific stuff
	//
	// Timestamp
	uint64_t a_timestamp_sec;
	uint64_t a_timestamp_nsec;

	// Frame properties
	int32_t a_x_size_max; // Sensor size
	int32_t a_y_size_max;
	int32_t a_x_start; // Selected size and positions
	int32_t a_x_end;
	int32_t a_y_start;
	int32_t a_y_end;
	int32_t a_x_bin;
	int32_t a_y_bin;

	// Camera settings
	int32_t a_ad_channel; // ADC
	int32_t a_amplifier; // EM or classical preamplifier
	int32_t a_preamp_gain; // Preamplifier gain
	int32_t a_em_gain; // EM gain
	int32_t a_hs_speed; // HS speed
	int32_t a_vs_speed; // VS speed
	int32_t a_vs_amplitude; // VS amplitude
	float a_exposure; // Exposure time in seconds
	int32_t a_shutter; // Shutter
	int32_t a_trigger; // Trigger
	int32_t a_temperature; // Temperature
	int32_t a_cooler; // Cooler
	int32_t a_cooler_mode; // Cooler mode
	int32_t a_fan; // Fan

	//
	// Image data
	//
	uint8_t *data;
} __attribute__((__packed__));

//
// Ximea movie frame struct
//
struct ximea_save_struct {
	//
	// Common stuff
	//
	uint32_t magic; // 'TemI'
	uint32_t version;
	uint32_t type; // Camera type
	uint32_t pixelmode; // Pixel mode
	uint32_t length_header; // Header data in bytes ( Everything except image data )
	uint32_t length_data; // Total data length in bytes;

	//
	// Camera specific stuff
	//
	// Timestamp
	uint64_t x_timestamp_sec;
	uint64_t x_timestamp_nsec;

	uint32_t x_size_x_max; // Sensor size
	uint32_t x_size_y_max;
	uint32_t x_size_x; // Selected region
	uint32_t x_size_y;
	uint32_t x_pos_x;
	uint32_t x_pos_y;

	// Camera settings

	//
	// Image data
	//
	uint8_t *data;
} __attribute__((__packed__));

#define NAMELENGTH 100

int main( int argc, char *argv[] )
{
	int i, j;
	FILE *moviefile;
	char moviename[NAMELENGTH], tiffname[NAMELENGTH];
	struct iidc_save_struct frame;
	uint8_t *imagebuf;
	TIFF *image;
	char *buf;
	int width, height;
	int start_index;
	double time_interval = 0.00033;  // fps  to add
	
	imagebuf = (uint8_t *)malloc( 2000*2000*2 );
	if ( imagebuf == NULL ) {
		printf( "Couldn't allocate enough memory.\n" );
		exit( EXIT_FAILURE );
	}
	
	if ( argc  != 3 ) {
		printf( "Wrong command.\n\"tiff2movie start_index basename\"\n" );
		exit( EXIT_FAILURE );
	}
	
	start_index = atoi( argv[1] );

	// Open a movie file
	sprintf( moviename, "%s.movie", argv[2] );
	if ( !( moviefile = fopen( moviename, "wb" ) ) ) {
		printf( "Couldn't create a movie file.\n" );
		exit( EXIT_FAILURE );
	}
	
	// Initialise the frame
	memset( &frame, 0, IIDC_MOVIE_HEADER_LENGTH );
	//
	// Common stuff
	//
	frame.magic = CAMERA_MOVIE_MAGIC; // 'TemI'
	frame.version = CAMERA_MOVIE_VERSION;
	frame.type = CAMERA_TYPE_IIDC; // Camera type
	frame.pixelmode = CAMERA_PIXELMODE_MONO_8 ; // Pixel mode
	frame.length_header = IIDC_MOVIE_HEADER_LENGTH; // Header data in bytes ( Everything except image data )
	//frame.length_data; // Total data length in bytes;
	
	buf = malloc( 20000 );
	
	// Read all the files
	i = start_index;
	sprintf( tiffname, "%s%d.tif", argv[2], i );
	while ( ( image = TIFFOpen( tiffname, "r" ) ) != NULL ) {
					
		// We need to set some values for basic tags before we can add any data
		/*TIFFSetField( image, TIFFTAG_IMAGEWIDTH, size_x );
		TIFFSetField( image, TIFFTAG_IMAGELENGTH, size_y );
		TIFFSetField( image, TIFFTAG_BITSPERSAMPLE, 8 * bpp );
		TIFFSetField( image, TIFFTAG_SAMPLESPERPIXEL, 1 );
		TIFFSetField( image, TIFFTAG_ROWSPERSTRIP, size_y );

		TIFFSetField( image, TIFFTAG_COMPRESSION, COMPRESSION_LZW );
		TIFFSetField( image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK );
		TIFFSetField( image, TIFFTAG_FILLORDER, FILLORDER_MSB2LSB );
		TIFFSetField( image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG );*/
		
		TIFFGetField( image, TIFFTAG_IMAGEWIDTH, &width );
		TIFFGetField( image, TIFFTAG_IMAGELENGTH, &height );
		frame.i_size_x = width;
		frame.i_size_y = height;
		frame.length_data = width * height;
		frame.i_data_depth = 8;
		frame.i_image_bytes = frame.length_data;
		frame.i_total_bytes = frame.length_data;
		frame.i_timestamp = time_interval * i *1e6;
		
		if ( fwrite( &frame, IIDC_MOVIE_HEADER_LENGTH, 1, moviefile ) != 1 ) {
			printf( "Error in writing a movie file\n" );
			exit( EXIT_FAILURE );
		}
		for ( j = 0; j < height; j++ ) {
			TIFFReadScanline( image , buf, j, 0 );
			if ( fwrite( buf, width, 1, moviefile ) != 1 ) {
				printf( "Error in writing a movie file\n" );
				exit( EXIT_FAILURE );
			}
		}
		// Close the file
		TIFFClose( image );	
				
		i++;
		sprintf( tiffname, "%s%d.tif", argv[2], i );
		if ( i % 1000 == 0 )
			printf( "frame %d\n", i );
	}
	
	return 0;
}
