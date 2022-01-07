#include "LibraryInfo.hpp"

#include <LibraryGitHash.hpp>
#include <LibraryVersion.hpp>

#include <string>

namespace PacBio {
namespace Actc {

Library::Bundle LibraryBundle()
{
    Library::Bundle bundle{LibraryInfo(), {}};
    return bundle;
}

Library::Info LibraryInfo() { return {"actc", Actc::ReleaseVersion, Actc::LibraryGitSha1}; }

}  // namespace Actc
}  // namespace PacBio
