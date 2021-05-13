using SpectrumData;
using System;
using System.Collections.Generic;
using System.Text;

namespace SpectrumProcess
{
    public interface IProcess
    {
        ISpectrum Process(ISpectrum spectrum);
    }
}
