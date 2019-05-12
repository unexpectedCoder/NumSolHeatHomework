#ifndef ERR_H
#define ERR_H

#include <string>


struct Error
{
  std::string mess = "";
  const std::string& sendEx(const std::string& ex)
  {
    mess = "\n\tError: " + ex + "!\n";
    return mess;
  }

  Error() {}
  ~Error()
  {
    mess.~basic_string();
  }
};


#endif // ERR_H
