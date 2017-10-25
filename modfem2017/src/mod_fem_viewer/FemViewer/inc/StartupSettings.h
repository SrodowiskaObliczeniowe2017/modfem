#ifndef _STARTUP_SETTINGS_H_
#define _STARTUP_SETTINGS_H_

#include "fv_config.h"

namespace FemViewer {

	// class for default sturtupo settings
	class StartupSettings {
	public:
		static const char* WorkDir;   // = work_dir;
		static const char* InDataDir; // = work_dir;
		static const char* OutDataDir;// = work_dir;
		static const char* FieldFile; // = field_file;
		static const char* MeshFile;//	  = mesh_file;
		static const char* FieldFileExt;// = file_ext;
		static const char* MeshFileExt;//	= file_ext;
		static const char* FieldFileMask;// = field_file_mask;
		static const char* MeshFileMask;//	 = mesh_file_mask;
	public:
		StartupSettings(){}
		~StartupSettings(){}
	};

	/*static*/ const char* StartupSettings::WorkDir    = work_dir;
	/*static*/ const char* StartupSettings::InDataDir  = work_dir;
	/*static*/ const char* StartupSettings::OutDataDir = work_dir;
	/*static*/ const char* StartupSettings::FieldFile  = field_file;
	/*static*/ const char* StartupSettings::MeshFile	  = mesh_file;
	/*static*/ const char* StartupSettings::FieldFileExt = file_ext;
	/*static*/ const char* StartupSettings::MeshFileExt	= file_ext;
	/*static*/ const char* StartupSettings::FieldFileMask = field_file_mask;
	/*static*/ const char* StartupSettings::MeshFileMask	 = mesh_file_mask;
} // end namespace FemViewer

#endif /* _STARTUP_SETTINGS_H_
*/