#include <csignal>
#include <functional>

namespace signal_handler {

/*
template <int SIGNAL>
struct signal_handler
{
  static std::function<void(int)> function_handler;

  static void
  signal ()
  {
    std::signal( SIGNAL ,
      [] (int signal) -> void {
        signal_handler<SIGNAL>::function_handler(SIGNAL);
        exit(-1);
      });
  }

  static int
  call ()
  { return std::raise(SIGNAL); }
};

template <int SIGNAL>
std::function<void(int)> signal_handler<SIGNAL>::function_handler;
*/

template <int SIGNAL>
struct signal_handler_one {
  static std::function<void(int)> function_handler;

  template <typename FUNC>
  static void
  handler ( FUNC f )
  {
    function_handler = f;
    // we use `std::signal` with one `int` (`SIGNAL`) and one `void(int)` function
    // this `void(int)` function can't capture any values, but `function_handler` can
    std::signal( SIGNAL ,
      [] (int signal) -> void {
        signal_handler_one<SIGNAL>::function_handler(signal);
        exit(-1);
    });
  }
};
template <int SIGNAL>
std::function<void(int)> signal_handler_one<SIGNAL>::function_handler;

template <int ...SIGNAL>
struct signal_handler {
  template < typename FUNC >
  static void
  handler ( FUNC function_handler )
  {
    int _[] = {0, (signal_handler_one<SIGNAL>::handler(function_handler),0)...}; // dummy array unpack at compilation time to loop over all signals in `...SIGNAL`
    (void)_; // beware of compilator optimisation
  }

  /**
    example:

    ```c++
      int i=42;
      signals_handler::signals_handler<SIGINT,SIGILL>::handler([&](int signal)->void {
        std::cerr << "what is the value of i ?\n...\n\n" << i << "\n";
      });
    ```
  **/

};

} // namespace signals_handler
