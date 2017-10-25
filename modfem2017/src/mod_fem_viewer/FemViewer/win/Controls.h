///////////////////////////////////////////////////////////////////////////////
// Controls.h
// ==========
// collection of controls (button, checkbox, textbox...)
//
// Button       : push button
// CheckBox     : check box
// RadioButton  : radio button
// TextBox      : static label box
// EditBox      : editable text box
// ListBox      : list box
// Trackbar     : trackbar(slider), required comctl32.dll
// ComboBox     : combo box
// TreeView     : tree-view control, required comctl32.dll
// UpDownBox    : Up-Down(Spin) control, required comctl32.dll
//
//  AUTHOR: Song Ho Ahn (song.ahn@gmail.com)
// CREATED: 2005-03-28
// UPDATED: 2008-09-15
///////////////////////////////////////////////////////////////////////////////

#ifndef WIN_CONTROLS_H
#define WIN_CONTROLS_H

#include <windows.h>
#include <commctrl.h>                   // common controls

namespace FemViewer {
namespace Win
{
    // constants //////////////////////////////////////////////////////////////
    enum { MAX_INDEX = 30000 };


    ///////////////////////////////////////////////////////////////////////////
    // base class of control object
    ///////////////////////////////////////////////////////////////////////////
    class ControlBase
    {
    public:
        // ctor / dtor
        ControlBase() : handle(0), parent(0), id(0), fontHandle(0) {}
        ControlBase(HWND parent, int id, bool visible=true) : handle(GetDlgItem(parent, id)) { if(!visible) disable(); }
        ~ControlBase() { if(fontHandle) DeleteObject(fontHandle);
                         fontHandle = 0; }

        // return handle
        HWND getHandle() const { return handle; }
        HWND getParentHandle() const { return parent; }

        // set all members
        void set(HWND parent, int id, bool visible=true) { this->parent = parent;
                                                           this->id = id;
                                                           handle = GetDlgItem(parent, id);
                                                           if(!visible) disable(); }

        // show/hide control
        void show() { ShowWindow(handle, SW_SHOW); }
        void hide() { ShowWindow(handle, SW_HIDE); }

        // set focus to the control
        void setFocus() { SetFocus(handle); }

        // enable/disable control
        void enable() { EnableWindow(handle, true); }
        void disable() { EnableWindow(handle, false); }
        bool isVisible() const { return (IsWindowVisible(handle) != 0); }

        // change font
        void setFont(char* fontName, int size, bool bold=false, bool italic=false, bool underline=false, unsigned char charSet=DEFAULT_CHARSET)
                  { HDC hdc = GetDC(handle);
                    LOGFONT logFont;
                    HFONT oldFont = (HFONT)SendMessage(handle, WM_GETFONT, 0 , 0);
                    GetObject(oldFont, sizeof(LOGFONT), &logFont);
                    strncpy(logFont.lfFaceName, fontName, LF_FACESIZE);
                    logFont.lfHeight = -MulDiv(size, GetDeviceCaps(hdc, LOGPIXELSY), 72);
                    ReleaseDC(handle, hdc);
                    if(bold) logFont.lfWeight = FW_BOLD;
                    else     logFont.lfWeight = FW_NORMAL;
                    logFont.lfItalic = italic;
                    logFont.lfUnderline = underline;
                    logFont.lfCharSet = charSet;

                    if(fontHandle) DeleteObject(fontHandle);
                    fontHandle = CreateFontIndirect(&logFont);
                    SendMessage(handle, WM_SETFONT, (WPARAM)fontHandle, MAKELPARAM(1, 0));
                  }


    protected:
        HWND  handle;
        HWND  parent;
        int   id;
        HFONT fontHandle;
    };



    ///////////////////////////////////////////////////////////////////////////
    // button class
    ///////////////////////////////////////////////////////////////////////////
    class Button : public ControlBase
    {
    public:
        // ctor / dtor
        Button() : ControlBase() {}
        Button(HWND parent, int id, bool visible=true) : ControlBase(parent, id, visible) {}
        ~Button() {}

        // change caption
        void setText(const wchar_t* text) { SendMessage(handle, WM_SETTEXT, 0, (LPARAM)text); }
    };



    ///////////////////////////////////////////////////////////////////////////
    // Checkbox class
    ///////////////////////////////////////////////////////////////////////////
    class CheckBox : public Button
    {
    public:
        CheckBox() : Button() {}
        CheckBox(HWND parent, int id, bool visible=true) : Button(parent, id, visible) {}
        ~CheckBox() {}

        void check() { SendMessage(handle, BM_SETCHECK, (WPARAM)BST_CHECKED, 0); }
        void uncheck() { SendMessage(handle, BM_SETCHECK, (WPARAM)BST_UNCHECKED, 0); }
        bool isChecked() const { return (SendMessage(handle, BM_GETCHECK, 0, 0) == BST_CHECKED); }
    };



    ///////////////////////////////////////////////////////////////////////////
    // Radio button class
    ///////////////////////////////////////////////////////////////////////////
    class RadioButton : public Button
    {
    public:
        RadioButton() : Button() {}
        RadioButton(HWND parent, int id, bool visible=true) : Button(parent, id, visible) {}
        ~RadioButton() {}

        void check() { SendMessage(handle, BM_SETCHECK, (WPARAM)BST_CHECKED, 0); }
        void uncheck() { SendMessage(handle, BM_SETCHECK, (WPARAM)BST_UNCHECKED, 0); }
        bool isChecked() const { return (SendMessage(handle, BM_GETCHECK, 0, 0) == BST_CHECKED); }
    };



    ///////////////////////////////////////////////////////////////////////////
    // Text box class
    ///////////////////////////////////////////////////////////////////////////
    class TextBox: public ControlBase
    {
    public:
        TextBox() : ControlBase() {}
        TextBox(HWND parent, int id, bool visible=true) : ControlBase(parent, id, visible) {}
        ~TextBox() {}

        void setText(const wchar_t* buf) { SendMessage(handle, WM_SETTEXT, 0, (LPARAM)buf); }
        void getText(wchar_t* buf, int len) const { SendMessage(handle, WM_GETTEXT, (WPARAM)len, (LPARAM)buf); }
        int  getTextLength() const { return (int)SendMessage(handle, WM_GETTEXTLENGTH, 0, 0); } // return number of chars
    };



    ///////////////////////////////////////////////////////////////////////////
    // Edit box class
    ///////////////////////////////////////////////////////////////////////////
    class EditBox: public TextBox
    {
    public:
        EditBox() : TextBox(), maxLength(0) {}
        EditBox(HWND parent, int id, bool visible=true) : TextBox(parent, id, visible), maxLength((int)SendMessage(handle, EM_GETLIMITTEXT, 0, 0)) {}
        ~EditBox() {}

        // set all members
        void set(HWND parent, int id, bool visible=true) { ControlBase::set(parent, id, visible);
                                                           maxLength = (int)SendMessage(handle, EM_GETLIMITTEXT, 0, 0); }

        void selectText() { SendMessage(handle, EM_SETSEL, 0, -1); }
        void unselectText() { SendMessage(handle, EM_SETSEL, -1, 0); }

        int getLimitText() const { return maxLength; } // return max number of chars

        static bool isChanged(int code) { return (code==EN_CHANGE); } // LOWORD(wParam)==id && HIWORD(wParam)==EN_CHANGE

    protected:
        int maxLength;
    };



    ///////////////////////////////////////////////////////////////////////////
    // List box class
    ///////////////////////////////////////////////////////////////////////////
    class ListBox: public ControlBase
    {
    public:
        ListBox() : ControlBase(), listCount(0) {}
        ListBox(HWND parent, int id, bool visible=true) : ControlBase(parent, id, visible), listCount(0) {}
        ~ListBox() {}

        int getCount() const { return listCount; }

        void resetContent() { SendMessage(handle, LB_RESETCONTENT, 0, 0); }

        // add a string at the end of the list
        // If list is full, then delete the most front entry first and add new string.
        // Also, make the latest entry is visible with LB_SETTOPINDEX.
        void addString(const wchar_t* str) { if(listCount >= MAX_INDEX) deleteString(0);
                                             int index = (int)SendMessage(handle, LB_ADDSTRING, 0, (LPARAM)str);
                                             if(index != LB_ERR) ++listCount;
                                             SendMessage(handle, LB_SETTOPINDEX, index, 0); }

        // insert a string at given index into the list
        // If the index is greater than listCount, the string is added to the end of the list.
        // If the list is full, then delete the last entry before inserting new string.
        void insertString(const wchar_t* str, int index) { if(listCount >= MAX_INDEX) deleteString(listCount-1);
                                                           index = (int)SendMessage(handle, LB_INSERTSTRING, index, (LPARAM)str);
                                                           if(index != LB_ERR) ++listCount;
                                                           SendMessage(handle, LB_SETTOPINDEX, index, 0); }

        void deleteString(int index) { if(SendMessage(handle, LB_DELETESTRING, index, 0) != LB_ERR) --listCount; }

    private:
        int listCount;
    };



    ///////////////////////////////////////////////////////////////////////////
    // Trackbar class (Slider)
    // It requires loading comctl32.dll. Note that the range values are logical
    // numbers and all integers.For example, if you want the range of 0~1 and
    // 100 increment steps(ticks) between 0 and 1, then you need to call
    // setRange(0, 100), not setRange(0,1). In other words, this class is not
    // responsible to scale the range to actual values.
    ///////////////////////////////////////////////////////////////////////////
    class Trackbar: public ControlBase
    {
    public:
        Trackbar() : ControlBase() {}
        Trackbar(HWND parent, int id, bool visible=true) : ControlBase(parent, id, visible) {}
        ~Trackbar() {}

        // set/get range
        void setRange(int first, int last) { SendMessage(handle, TBM_SETRANGE, (WPARAM)true, (LPARAM)MAKELONG(first, last)); }
        int getRangeMin() const { return (int)SendMessage(handle, TBM_GETRANGEMIN, 0, 0); }
        int getRangeMax() const { return (int)SendMessage(handle, TBM_GETRANGEMAX, 0, 0); }

        // set a tick mark at a specific position
        // Trackbar creates its own first and last tick marks, so do not use it for first and last tick mark.
        void setTic(int pos) { if(pos <= getRangeMin()) return;   // skip if it is the first tick
                               if(pos >= getRangeMax()) return;   // skip if it is the last tick
                               SendMessage(handle, TBM_SETTIC, 0, (LPARAM)pos); }

        // set the interval frequency for tick marks.
        void setTicFreq(int freq) { SendMessage(handle, TBM_SETTICFREQ, (WPARAM)freq, 0); }

        // set/get current thumb position
        int  getPos() const { SendMessage(handle, TBM_GETPOS, 0, 0); }
        void setPos(int pos) { SendMessage(handle, TBM_SETPOS, (WPARAM)true, (LPARAM)pos); }
    };



    ///////////////////////////////////////////////////////////////////////////
    // Combo box class
    ///////////////////////////////////////////////////////////////////////////
    class ComboBox: public ControlBase
    {
    public:
        ComboBox() : ControlBase() {}
        ComboBox(HWND parent, int id, bool visible=true) : ControlBase(parent, id, visible) {}
        ~ComboBox() {}

        int getCount() const { return (int)SendMessage(handle, CB_GETCOUNT, 0, 0); }

        void resetContent() { SendMessage(handle, CB_RESETCONTENT, 0, 0); }

        // add a string at the end of the combo box
        void addString(const wchar_t* str) { SendMessage(handle, CB_ADDSTRING, 0, (LPARAM)str); }

        // insert a string at given index into the list
        void insertString(const wchar_t* str, int index) { SendMessage(handle, CB_INSERTSTRING, index, (LPARAM)str); }

        // delete an entry
        void deleteString(int index) { SendMessage(handle, CB_DELETESTRING, index, 0); }

        // get/set current selection
        int getCurrentSelection() { return (int)SendMessage(handle, CB_GETCURSEL, 0, 0); }
        void setCurrentSelection(int index) { SendMessage(handle, CB_SETCURSEL, index, 0); }
    };






    ///////////////////////////////////////////////////////////////////////////
    // UpDownBox class
    // It requires comctl32.dll.
    // Note that UpDownBox sends UDN_DELTAPOS notification when the position of
    // the control is about to change. Return non-zero value to prevent the
    // change.
    ///////////////////////////////////////////////////////////////////////////
    class UpDownBox: public ControlBase
    {
    public:
        UpDownBox() : ControlBase() {}
        UpDownBox(HWND parent, int id, bool visible=true) : ControlBase(parent, id, visible) {}
        ~UpDownBox() {}

        // set/get the buddy window for up-down control
        // UDS_SETBUDDYINT style must be added to make automatic update of buddy window.
        void setBuddy(HWND buddy)                   { SendMessage(handle, UDM_SETBUDDY, (WPARAM)buddy, 0); }
        HWND getBuddy()                             { return (HWND)SendMessage(handle, UDM_GETBUDDY, 0, 0); }

        // set/get the base(decimal or hexadecimal) of number in the buddy window
        // 10 for deciaml, 16 for hexadecimal
        void setBase(int base)                      { SendMessage(handle, UDM_SETBASE, (WPARAM)base, 0); }
        int  getBase()                              { return (int)SendMessage(handle, UDM_GETBASE, 0, 0); }

        // set/get the range (low and high value) of up-down control
        void setRange(int low, int high)            { SendMessage(handle, UDM_SETRANGE32, (WPARAM)low, (LPARAM)high); }
        void getRange(int* low, int* high)          { SendMessage(handle, UDM_GETRANGE32, (WPARAM)low, (LPARAM)high); }

        // set/get acceleration (elapsed time and increment) of up-down control
        void setAccel(unsigned int second, unsigned int increment) { UDACCEL accels[2];
                                                                     accels[0].nSec = 0;      accels[0].nInc = 1;
                                                                     accels[1].nSec = second; accels[1].nInc = increment;
                                                                     SendMessage(handle, UDM_SETACCEL, (WPARAM)2, (LPARAM)&accels); }
        void getAccel(unsigned int* second, unsigned int* increment) { UDACCEL accels[2];
                                                                       SendMessage(handle, UDM_GETACCEL, (WPARAM)2, (LPARAM)&accels);
                                                                       *second = accels[1].nSec; *increment = accels[1].nInc; }

        // set/get the current position with 32-bit precision
        void setPos(int position)                   { SendMessage(handle, UDM_SETPOS32, 0, (LPARAM)position); }
        int  getPos()                               { return (int)SendMessage(handle, UDM_GETPOS32, 0, 0); }
    };

}
}

#endif