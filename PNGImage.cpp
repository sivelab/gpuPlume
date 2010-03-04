#include <cmath>
#include <cstdlib>
#include "PNGImage.h"

#define PROGNAME "pbuffer"
#define TOO_SMALL 0
#define BIG_ENOUGH 1

using namespace cs5721;

bool PNGImage::writeFileData(const std::string& filename, const int width, const int height, const float *data)
{
  int rc, error = 0;
  ulg rowbytes;

  m_png_fileinfo.image_data = NULL;
  m_png_fileinfo.row_pointers = NULL;
  m_png_fileinfo.filter = FALSE;
  m_png_fileinfo.interlaced = FALSE;
  m_png_fileinfo.have_bg = FALSE;
  m_png_fileinfo.have_time = FALSE;
  m_png_fileinfo.have_text = 0;
  m_png_fileinfo.gamma = 0;
  m_png_fileinfo.width = width;
  m_png_fileinfo.height = height;
  m_png_fileinfo.sample_depth = 8;
  m_png_fileinfo.modtime = time(0);
  m_png_fileinfo.infile = NULL;

  m_png_fileinfo.pnmtype = 6;

  /* check if outname already exists; if not, open */
  /* if ((m_png_fileinfo.outfile = fopen(outname, "rb")) != NULL) {
     fprintf(stderr, "  output file exists [%s]\n", outname);
     fclose(m_png_fileinfo.outfile);
     } else */
  if (!(m_png_fileinfo.outfile = fopen(filename.c_str(), "wb"))) {
    fprintf(stderr, "  can't open output file [%s]\n", filename.c_str());
  }
  
  /* allocate libpng stuff, initialize transformations, write pre-IDAT data */
  if ((rc = writepng_init(&m_png_fileinfo)) != 0) {
    switch (rc) {
      case 2:
	fprintf(stderr, PROGNAME
		":  libpng initialization problem (longjmp)\n");
	break;
      case 4:
	fprintf(stderr, PROGNAME ":  insufficient memory\n");
	break;
      case 11:
	fprintf(stderr, PROGNAME
		":  internal logic error (unexpected PNM type)\n");
	break;
      default:
	fprintf(stderr, PROGNAME
		":  unknown writepng_init() error\n");
	break;
    }
    exit(rc);
  }

  /* calculate rowbytes on basis of image type; note that this becomes
   * much more complicated if we choose to support PBM type, ASCII PNM
   * types, or 16-bit-per-sample binary data [currently not an
   * official NetPBM type] */

  if (m_png_fileinfo.pnmtype == 5)
    rowbytes = m_png_fileinfo.width;
  else if (m_png_fileinfo.pnmtype == 6)
    rowbytes = m_png_fileinfo.width * 3;
  else /* if (m_png_fileinfo.pnmtype == 8) */
    rowbytes = m_png_fileinfo.width * 4;

  /* read and write the image, either in its entirety (if writing
   * interlaced PNG) or row by row (if non-interlaced) */

  // fprintf(stderr, "Encoding image data...\n");
  // fflush(stderr);
  
  if (m_png_fileinfo.interlaced) {
    long i;
    ulg bytes;
    ulg image_bytes = rowbytes * m_png_fileinfo.height;   /* overflow? */

    m_png_fileinfo.image_data = (uch *)malloc(image_bytes);
    m_png_fileinfo.row_pointers = (uch **)malloc(m_png_fileinfo.height*sizeof(uch *));
    if (m_png_fileinfo.image_data == NULL || m_png_fileinfo.row_pointers == NULL) {
      fprintf(stderr, PROGNAME ":  insufficient memory for image data\n");
      writepng_cleanup(&m_png_fileinfo);
      cleanup();
      exit(5);
    }
    for (i = 0;  i < m_png_fileinfo.height;  ++i)
      m_png_fileinfo.row_pointers[i] = m_png_fileinfo.image_data + i*rowbytes;
    bytes = fread(m_png_fileinfo.image_data, 1, image_bytes, m_png_fileinfo.infile);
    if (bytes != image_bytes) {
      fprintf(stderr, PROGNAME ":  expected %lu bytes, got %lu bytes\n",
	      image_bytes, bytes);
      fprintf(stderr, "  (continuing anyway)\n");
    }
    if (writepng_encode_image(&m_png_fileinfo) != 0) {
      fprintf(stderr, PROGNAME
	      ":  libpng problem (longjmp) while writing image data\n");
      writepng_cleanup(&m_png_fileinfo);
      cleanup();
      exit(2);
    }
    
  } else /* not interlaced:  write progressively (row by row) */ {
    long j,x;
    ulg bytes;

    m_png_fileinfo.image_data = (uch *)malloc(rowbytes);
    if (m_png_fileinfo.image_data == NULL) {
      fprintf(stderr, PROGNAME ":  insufficient memory for row data\n");
      writepng_cleanup(&m_png_fileinfo);
      cleanup();
      exit(5);
    }
    error = 0;
    for (j = m_png_fileinfo.height-1;  j >= 0L;  --j) {
      long pixelcount;
      bytes = rowbytes;
      pixelcount = 0;
      for (x=0; x<rowbytes; x+=3) {
	m_png_fileinfo.image_data[x] = (uch)(floor(data[j*width*3 + pixelcount + 0] * 255));
	m_png_fileinfo.image_data[x+1] = (uch)(floor(data[j*width*3 + pixelcount + 1] * 255));
	m_png_fileinfo.image_data[x+2] = (uch)(floor(data[j*width*3 + pixelcount + 2] * 255));
	pixelcount+=3;
	/* printf("rowbytes=%ld, x=%ld, color = [%d, %d, %d]\n", rowbytes, x, m_png_fileinfo.image_data[x+0], m_png_fileinfo.image_data[x+1], m_png_fileinfo.image_data[x+2]); */
      }

      if (bytes != rowbytes) {
	fprintf(stderr, PROGNAME
		":  expected %lu bytes, got %lu bytes (row %ld)\n", rowbytes,
		bytes, m_png_fileinfo.height-j);
	++error;
	break;
      }
      if (writepng_encode_row(&m_png_fileinfo) != 0) {
	fprintf(stderr, PROGNAME
		":  libpng problem (longjmp) while writing row %ld\n",
		m_png_fileinfo.height-j);
	++error;
	break;
      }
    }
    if (error) {
      writepng_cleanup(&m_png_fileinfo);
      cleanup();
      exit(2);
    }
    if (writepng_encode_finish(&m_png_fileinfo) != 0) {
      fprintf(stderr, PROGNAME ":  error on final libpng call\n");
      writepng_cleanup(&m_png_fileinfo);
      cleanup();
      exit(2);
    }
  }
    
  /* OK, we're done (successfully):  clean up all resources and quit */

  // fprintf(stderr, "Done.\n");
  // fflush(stderr);

  writepng_cleanup(&m_png_fileinfo);
  return true;
}


void PNGImage::cleanup()
{
  if (m_png_fileinfo.outfile) {
    fclose(m_png_fileinfo.outfile);
    m_png_fileinfo.outfile = NULL;
  }

  if (m_png_fileinfo.infile) {
    fclose(m_png_fileinfo.infile);
    m_png_fileinfo.infile = NULL;
  }

  if (m_png_fileinfo.image_data) {
    free(m_png_fileinfo.image_data);
    m_png_fileinfo.image_data = NULL;
  }

  if (m_png_fileinfo.row_pointers) {
    free(m_png_fileinfo.row_pointers);
    m_png_fileinfo.row_pointers = NULL;
  }
}
