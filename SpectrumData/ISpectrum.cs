using System.Collections.Generic;

namespace SpectrumData
{
    public enum TypeOfMSActivation
    { CID, MPD, ECD, PQD, ETD, HCD, Any, SA, PTR, NETD, NPTR };
    public interface ISpectrum
    {
        List<IPeak> GetPeaks();
        void SetPeaks(List<IPeak> peaks);
        void Add(IPeak peak);
        void Clear();
        int GetScanNum();
        double GetRetention();
        void set_scan(int scan);
        void set_retention(double retention);

        ISpectrum Clone();
    }
}
