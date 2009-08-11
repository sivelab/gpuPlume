#ifndef RTC_SYNC_H
#define RTC_SYNC_H

#include <linux/rtc.h>
#include <sys/ioctl.h>
#include <sys/time.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <stdio.h>

class rtcSync 
{
public:
  rtcSync( unsigned int hz = 32 );
  ~rtcSync();

  void block();

private:
  int _dev_fd;
  unsigned int _hz;
};

inline void rtcSync::block()
{
  unsigned long data;  
  read( _dev_fd, &data, sizeof(unsigned long) );
}

#endif // RTC_SYNC_H
