import numpy as np
import subprocess
import glob
import os

# dspsr -A -b 128 -L 10 -scloffs -c "0.08990380293086397554" -O J1326-4728 -e cf -U 900 -t 10 uwl_220213_144049.sf
all_cal = glob.glob("*.cf")

for filename in all_cal:
    cmd = 'psredit -m -c type="PolnCal" %s'%filename
    print (cmd)
    subprocess.call(cmd, shell=True)

    cmd = 'psredit -c rcvr:basis=lin,rcvr:hand=-1,rcvr:sa=0.0,rcvr:rph=0.0 -m %s'%filename
    print (cmd)
    subprocess.call(cmd, shell=True)

    print (filename)
    output = filename + '.dzT'
    zapped = filename + '.pazi'

    if os.path.isfile(output):
        print ('%s done!'%output)
    else:
        cmd = 'paz -F "749 797" -F "863 898" -F "918 932" -F "938 968" -F "1014 1031" -F "1077 1092" -F "1120 1122" -F "1159 1178" -F "1260 1290" -F "1553 1557" -F "1557 1562" -F "1574 1577" -F "1581 1582" -F "1620 1627" -F "1728 1734" -F "1741 1746" -F "1805 1824" -F "1845 1867" -F "1919 1929" -F "1975 1984" -F "2140 2150" -F "2328 2362" -F "2409 2477" -F "3447 3457" -F "3552 3569" -F "699 709" -F "827 837" -F "955 965" -F "1083 1093" -F "1211 1221" -F "1339 1349" -F "1467 1477" -F "1595 1605" -F "1723 1733" -F "1851 1861" -F "1979 1989" -F "2107 2117" -F "2235 2245" -F "2363 2373" -F "2491 2501" -F "2619 2629" -F "2747 2757" -F "2875 2885" -F "3003 3013" -F "3131 3141" -F "3259 3269" -F "3387 3397" -F "3515 3525" -F "3643 3653" -F "3771 3781" -F "3899 3909" -F "4027 4037" -z "1926 1927 1928 1929 1930 1931 1932 1933 1934 1935 1936 1937 1938 1939 1940 1941 1942 1943 1944 1122 1123 1124 1125 1126 1127 1128 1129 1130 1131 1132 1133 1134 1135 1136 1137 1138 1139 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 5 6 7 8 9 10 28 29 213 100 132 133 134 135 136 137 138 139 140 211 212 1618 1619 1620 1621 1622 1623 1624 1625 1626 1627 1628 1629 1630 1631 1632 1633 1634 1635 1636 1637 1638 1639 1640 1641 1642 1643 1644 1645 1646 1647 1648 1649 1650 1651 1652 1653 1654 1655 1656 1657 1658 1659 1660 1661 1662 1663 1664 1665 1666 1667 1668 1669 1670 1671 1672 1673 1674 1675 1676 1773 1774 1775 1776 1777 1697 1698 1699 1832 1833 1834 1835 1836 1837 1838 1839 1840 1841 1842 1843 1844 1845 1701 1702 1703 1704 1705 1706 1707 337 336 327 328 388 389 390 391 392 393 394 504 505 506 507 422 423 424 425 426 427 428 429 430 431 432 433 434 120 121 122 123 1451 1452 1453 1454 1455 1456 1457 1458 1459 1460 1461 1462 1463 1943 1944 1945 1946 1947 1948 1949 1950 1951 1952 1953 1954 1955 1956 1957 1958 1959 1960 1961 1962 1963 1964 1965 1966 1967 1968 1969 1970 1971 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985 1986 2871 2872 2873 2874 2875 2876 2877 2878 2879 2880 2881 2882 2883 2884 2885 2886 2887 2888 2889 2890 2891 2892 2893 2894 2895 2896 2897 2898 2899 2900 2901 2902 2903 2904 2905 2906 2907 2908 2909 2910 2911 2912 2913 2914 2915 2916 2917 2918 2919 2920 2921 2922 2923 2924 2925 2926 2927 2928 2929 2930 2740 2741 2742 2743 2744 2745 2746 2747 2748 2749 2750 2751 2752 2753 2754 2755 2756 2757 2758 2759 2760 2761 2846 2847 2848 2849 2850 914 915 1119 1120 1121 1122 1123 1029 1030 1031 1032 1006 1007 1008 1009 1010 329 330 331 332" -w "0 12" -e cf.pazi %s'%filename
        print (cmd)
        subprocess.call(cmd, shell=True)

        cmd = 'pam -T -e dzT %s'%zapped
        print (cmd)
        subprocess.call(cmd, shell=True)

        print ('Cleaning %s'%filename)