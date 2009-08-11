#include "rtcSync.h"

rtcSync::rtcSync( unsigned int hz )
  : _hz( hz )
{
  char *rtc_devname = (char *)"/dev/rtc";

  if( (_dev_fd = open( rtc_devname, O_RDONLY )) < 0 )
    {
      char emsg[128];
      sprintf( emsg, "RTfs Device \"%s\"", rtc_devname );
      perror( emsg );
    }

  if( ioctl( _dev_fd, RTC_IRQP_SET, _hz ) < 0 )
    {
      char emsg[128];
      sprintf( emsg, "set update rate(%d)", _hz );
      perror( emsg );
      switch( errno )
        {
	  case EINVAL :
	    fprintf( stderr,
		     "(Hz) must be a power of 2 between 2 and 8192, inclusive.\n" );
	    break;

	  case EBADF :
	    fprintf( stderr,
                    "Device probably not initialized successfully.\n"
                    "  1. Check permissions on %s\n"
                    "  2. See if other processes may have control of %s\n"
                    "  3. Check environmental variable \"RTCDEV\"\n"
		     , rtc_devname, rtc_devname );
	    break;

	  case EACCES :
	    fprintf( stderr,
		     "You must have root priviledges to run at update rates higher than 64 hz.\n" );
	    break;

	  default :
	    break;
        }
    }

  if( ioctl( _dev_fd, RTC_PIE_ON, 0 ) < 0 )
    perror( "cannot start periodic interrupts" );
}

rtcSync::~rtcSync()
{
  int retval = ioctl( _dev_fd, RTC_PIE_OFF, 0 );
  if (retval == -1) {
    perror("ioctl");
  }

  close( _dev_fd );
}
