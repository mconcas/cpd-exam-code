#include "util.hpp"
#include "cl.hpp"

using std::cout;
using std::vector;
using std::string;
using std::endl;
using std::size_t;

inline string loadProgram(string input);
void show_platforms(vector<cl::Platform> &platforms);

int main() {
  int err;
  std::vector<cl::Platform> platforms;
  err = cl::Platform::get(&platforms);
  if (err == CL_INVALID_VALUE) return -666;
  show_platforms(platforms);
  cl::Context contesto(CL_DEVICE_TYPE_GPU);
  // cl::Program programma(contesto, loadProgram("routines/Square.cl"));
  auto devices = contesto.getInfo<CL_CONTEXT_DEVICES>();
  // programma.build(devices);
  return err;
}

inline std::string loadProgram(std::string input)
{
  std::ifstream stream(input.c_str());
  if (!stream.is_open()) {
    std::cout << "Cannot open file: " << input << std::endl;
    exit(1);
  }
  return std::string(std::istreambuf_iterator<char>(stream),
    (std::istreambuf_iterator<char>()));
}

void show_platforms(vector<cl::Platform> &platforms)
{
  for(auto platform: platforms) {
    vector<cl_platform_info > names = { CL_PLATFORM_EXTENSIONS,
			       CL_PLATFORM_NAME,
                               CL_PLATFORM_PROFILE,
                               CL_PLATFORM_VENDOR,
                               CL_PLATFORM_VERSION };
    vector<string> keys = {"extensions", "name", "profile", "vendor", "version"};
    cl::STRING_CLASS param;

    cout<<"Platform informations:"<<endl;
    for(size_t i=0; i<names.size(); ++i) {
      platform.cl::Platform::getInfo(names[i], &param);
      // auto pos = std::distance(names.begin(), name);
      cout<<"\t "<< keys[i] <<": "<< param <<endl;
    }
  }
}
