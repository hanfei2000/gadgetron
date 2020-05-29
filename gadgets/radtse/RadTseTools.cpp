#include "RadTseTools.h"

void Gadgetron::RadTse::newViewSelect(std::vector<std::vector<int>> &view, int nviews, int etl) {

    int num_et, number;
    int i, e;

    //calculate 2 sets of pseudo GA orders, one for lNumEchoTrains (or TRs), another for lEchoTrainLength (or TEs)
    long alPGA_tr[1000], alPGA_adc[1000], alDupLoc[1000], alMisNum[1000]; //arrays (vector) needed, with plenty of space
    double dPIcover = 1.0;
    if ((nviews % 2) == 1)
        dPIcover = 2.0;
    double dGA = 111.246117975 / (dPIcover * 180.0); //golden angle ratios
    long k = 0, m = 0, dup = 0, mis = 0;
    bool bFind = false;

    num_et = (int)std::ceil((double)nviews / etl);

    //first, calculate alPGA_tr (lNumEchoTrains=num_et)
    for (k = 0; k < num_et; k++)
    {
        alPGA_tr[k] = (long)((k * dGA) * num_et) % num_et;
        if (alPGA_tr[k] == num_et)
            alPGA_tr[k] = 0;
    }
    //find duplicates and replace it with missing numbers
    for (k = 0; k < num_et; k++)
    {
        bFind = false; //reset bfind
        for (m = 0; m < num_et; m++)
        {
            if (k == alPGA_tr[m])
            {
                if (!bFind)
                {
                    bFind = true;
                }
                else
                {
                    alDupLoc[dup] = m;
                    dup += 1;
                }
            }
        }
        if (!bFind)
        {
            alMisNum[mis] = k;
            mis += 1;
        }
    }
    for (k = 0; k < dup; k++)
    {
        alPGA_tr[alDupLoc[k]] = alMisNum[k];
    }

    //repeat this for alPGA_adc (m_lEchoTrainLength=etl)
    for (k = 0; k < etl; k++)
    {
        alPGA_adc[k] = (long)((k * dGA) * etl) % etl;
        if (alPGA_adc[k] == etl)
            alPGA_adc[k] = 0;
    }
    //find duplicates and replace with missing numbers
    dup = 0;
    mis = 0;
    for (k = 0; k < etl; k++)
    {
        bFind = false; //reset bfind
        for (m = 0; m < etl; m++)
        {
            if (k == alPGA_adc[m])
            {
                if (!bFind)
                {
                    bFind = true;
                }
                else
                {
                    alDupLoc[dup] = m;
                    dup += 1;
                }
            }
        }
        if (!bFind)
        {
            alMisNum[mis] = k;
            mis += 1;
        }
    }
    for (k = 0; k < dup; k++)
    {
        alPGA_adc[alDupLoc[k]] = alMisNum[k];
    }
    number = (long)(std::floor((float)num_et / (float)etl + 0.5)) * etl; //number = N/ETL rounded up to the nearest multiple of ETL
    if (number < etl)
    {
        number = etl;
    } //minimum is ETL
    for (e = 0; e < num_et; e++)
    {
        for (i = 0; i < etl; i++)
        {
            view[i][alPGA_tr[e]] = ((i * number + alPGA_adc[i] + alPGA_tr[e] * (etl)) % nviews);
        }
    }

}


long Gadgetron::RadTse::prepareData(double *mask, int ncol, int nviews, int etl, int echoNum, double *denscomp, double USFactor)
{
    int tier, num_et;
    int i, echoIdx, echoTrainIdx, dummyIdx;
    int ctr, found, count(0);

    num_et = nviews / etl;
    std::vector<int> totEchos(etl);
    std::vector<int> rho(etl);

    for (i = 0; i < etl - 1; i++)
        rho[i] = 10 * i;
    rho[etl - 1] = ncol;

    std::vector<std::vector<int>> keepList, view;
    keepList.resize(etl);
    for (i = 0; i < etl; i++)
        keepList[i].resize(etl);

    view.resize(etl);
    for (i = 0; i < etl; i++)
        view[i].resize(num_et);

    Gadgetron::RadTse::newViewSelect(view, nviews, etl);

    /* Generate a list of views to keep for each tier */
    keepList[0][0] = echoNum;
    totEchos[0] = 1;

    for (tier = 1; tier < etl - 1; tier++)
    {
        /* All tiers should have this echo */
        totEchos[tier] = 0;
        keepList[tier][0] = echoNum;
        totEchos[tier]++;

        ctr = 1;
        while (ctr < (tier + 1))
        {
            if ((echoNum - ctr >= 0) && (echoNum - ctr < etl))
            {
                keepList[tier][totEchos[tier]] = (echoNum - ctr);
                totEchos[tier]++;
            }
            if ((echoNum + ctr >= 0) && (echoNum + ctr < etl))
            {
                keepList[tier][totEchos[tier]] = (echoNum + ctr);
                totEchos[tier]++;
            }
            ctr++;
        }
    }

    /* Last tier has all the echos */
    for (echoIdx = 0; echoIdx < etl; echoIdx++)
    {
        keepList[etl - 1][echoIdx] = echoIdx;
    }

    for (tier = 0; tier < etl - 1; tier++)
    {
        rho[tier] = (int)std::floor((double)((nviews / etl) / M_PI * USFactor) * totEchos[tier]); //z003cnsa...prototype

        if (rho[tier] > ncol / 2)
        {
            /* Rho can not be larger than the available data. Truncate. */
            rho[tier] = ncol / 2;
        }
        /*printf(" %d ",rho[tier]);*/
    }
    /*printf("\n totEchos[7]=%d",totEchos[7]);*/

    for (echoTrainIdx = 0; echoTrainIdx < num_et; echoTrainIdx++)
    {
        for (echoIdx = 0; echoIdx < etl; echoIdx++)
        {

            for (tier = 0; tier < etl - 1; tier++)
            {

                /* Delete if not in keepList */
                found = 0;
                for (dummyIdx = 0; dummyIdx < totEchos[tier]; dummyIdx++)
                {
                    if (keepList[tier][dummyIdx] == echoIdx)
                    {
                        found = 1;
                    }
                }

                switch (tier)
                {
                case 0:

                    if (!found)
                    {
                        for (i = ncol / 2 - rho[tier]; i <= ncol / 2 + rho[tier]; i++)
                        {
                            mask[view[echoIdx][echoTrainIdx] * ncol + i] = 1.0;
                            count++;
                        }
                    }
                    break;
                default:
                    if (!found)
                    {

                        if (rho[tier] == ncol / 2)
                        {
                            /* Nyquist radius larger than available	data. Use all data. */
                            /* If the previous tier was all above the max radius. No need to delete */
                            /* the same points again. */
                            if (rho[tier - 1] != ncol / 2)
                            {
                                /* The previous	tier was not above the max radius. */

                                /* Left	side */
                                for (i = 0; i <= ncol / 2 - rho[tier - 1] - 1; i++)
                                {
                                    mask[view[echoIdx][echoTrainIdx] * ncol + i] = 1.0;
                                    count++;
                                }
                                /* Right side */
                                for (i = ncol / 2 + rho[tier - 1] + 1; i < ncol; i++)
                                {
                                    mask[view[echoIdx][echoTrainIdx] * ncol + i] = 1.0;
                                    count++;
                                }
                            }
                        }
                        else
                        {

                            /* Left	side */
                            for (i = ncol / 2 - rho[tier]; i <= ncol / 2 - rho[tier - 1] - 1; i++)
                            {
                                mask[view[echoIdx][echoTrainIdx] * ncol + i] = 1.0;
                                count++;
                            }
                            /* Right side */
                            for (i = ncol / 2 + rho[tier - 1] + 1; i <= ncol / 2 + rho[tier]; i++)
                            {
                                mask[view[echoIdx][echoTrainIdx] * ncol + i] = 1.0;
                                count++;
                            }
                        }
                    }
                    break;
                }
            }
        }
    }

    if (denscomp != NULL)
    {
        int shot = nviews / etl;
        int VIEW;
        for (int i = 0; i < ncol; i++)
        {
            int t = 0;
            while (rho[t] < abs(i - ncol / 2))
            {
                t++;
            }
            if (totEchos[t] == 0)
                VIEW = nviews;
            else
                VIEW = totEchos[t] * shot;

            for (int j = 0; j < nviews; j++)
            {
                if (i == ncol / 2)
                {
                    denscomp[j * ncol + i] = M_PI / (4 * VIEW);
                }
                else
                {
                    denscomp[j * ncol + i] = M_PI * abs(i - ncol / 2) / VIEW;
                }
            }
        }
    }
    return count;
}
