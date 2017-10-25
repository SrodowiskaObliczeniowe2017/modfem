/*
 * GraphicElem.hpp
 *
 *  Created on: 2011-03-20
 *      Author: pawel
 */
#include "GraphicElem.hpp"


namespace FemViewer
{
    HostTransferFunction::HostTransferFunction() :
        m_breakpoints(),
        m_localDeviceTransferFunction(),
        m_dirty(true)
    {
        m_localDeviceTransferFunction.Init();
    }


    void HostTransferFunction::FreeDeviceMemory()
    {

    }

    void HostTransferFunction::AllocateDeviceMemory()
    {

    }
/*
    void HostTransferFunction::CopyToOptix(optixu::Context context, OptiXBuffer<ElVisFloat>& breakpoints, OptiXBuffer<ElVisFloat>& values, TransferFunctionChannel channel)
    {
        breakpoints.SetContext(context);
        breakpoints.SetDimensions(m_breakpoints.size());
        values.SetContext(context);
        values.SetDimensions(m_breakpoints.size());
        BOOST_AUTO(breakpointData, breakpoints.Map());
        BOOST_AUTO(valueData,values.Map());

        int index = 0;
        for(std::map<double, Breakpoint>::iterator iter = m_breakpoints.begin(); iter != m_breakpoints.end(); ++iter)
        {
            breakpointData[index] = (*iter).first;

            if( channel == eDensity ) valueData[index] = (*iter).second.Density;
            if( channel == eRed ) valueData[index] = (*iter).second.Col.Red();
            if( channel == eGreen ) valueData[index] = (*iter).second.Col.Green();
            if( channel == eBlue ) valueData[index] = (*iter).second.Col.Blue();

            ++index;
        }
    }
*/

    void HostTransferFunction::CopyToDeviceMemory()
    {

    }

    void HostTransferFunction::SynchronizeDeviceIfNeeded()
    {
        if (!m_dirty) return;

        FreeDeviceMemory();
        AllocateDeviceMemory();
        CopyToDeviceMemory();

        // Reset after OptiX updates.
        m_dirty = true;
    }

    void HostTransferFunction::SynchronizeOptiXIfNeeded()
    {
        if( !m_dirty) return;
    }

    void HostTransferFunction::SetBreakpoint(double s, const ColorRGBA& c)
    {
        SetBreakpoint(s, c, c.A);
    }

    void HostTransferFunction::SetBreakpoint(double s, const ColorRGBA& c, const mfvFloat_t& density)
    {
        if( m_breakpoints.find(s) == m_breakpoints.end() )
        {
            m_breakpoints[s].Col = c;
            m_breakpoints[s].Density = density;
            m_breakpoints[s].Scalar = static_cast<mfvFloat_t>(s);
            m_dirty = true;
            //OnTransferFunctionChanged();
        }
    }

    void HostTransferFunction::Clear()
    {
        if( m_breakpoints.size() > 0 )
        {
            m_breakpoints = std::map<double, Breakpoint>();
            m_dirty = true;
            //OnTransferFunctionChanged();
        }
    }

} // end namespace FemViewer
