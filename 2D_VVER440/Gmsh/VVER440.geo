// Gmsh project created on Thu Jun 27 10:28:22 2013
lc = 100 ;
cm =1 ;    // 0.01 
//   cell  1 
l =8.487048957 * cm ;
ll = 25.46114687* cm ; // ll = l*3
r =7.35* cm ;
rr =14.7   * cm ;  // rr = r*2
rr1 =14.7   * cm ;  // rr = r*2
rr2 =14.7   * cm ;  // rr = r*2
rr3 =14.7   * cm ;  // rr = r*2
rr4 =14.7   * cm ;  // rr = r*2
rr5 =14.7   * cm ;  // rr = r*2
rr6 =14.7   * cm ;  // rr = r*2
rr7 =14.7   * cm ;  // rr = r*2
rr8 =14.7   * cm ;  // rr = r*2
rr9 =14.7   * cm ;  // rr = r*2
rr10 =14.7   * cm ;  // rr = r*2
rr11 =14.7   * cm ;  // rr = r*2
rr12 =14.7   * cm ;  // rr = r*2
rr13 =14.7   * cm ;  // rr = r*2
rr14 =14.7   * cm ;  // rr = r*2
rr15 =14.7   * cm ;  // rr = r*2


rm =14.7   * cm ;  // rr = r*2
  rmm = 14.7  ;
  rm2 = 14.7  ;
  rm3 = 14.7  ;
  rm4 = 14.7  ;
  rm5 = 14.7  ;
  rm6 = 14.7  ;
  rm7 = 14.7  ;
  rm8 = 14.7  ;
  rm9 = 14.7  ;
  rm10 = 14.7  ;
  rm11 = 14.7  ;
  rm12 = 14.7  ;
  rm13 = 14.7  ;
  rm14 = 14.7  ;
  rm15 = 14.7  ;
s = 4.243524478* cm ;  // R- 0.5 l = 0.5 *l
ls = 12.730573435 * cm ;   // ll = l+s
Point(1) = {0, 0, 0, lc} ;
Point(2) = {r, -s, 0, lc} ;
Point(3) = {rr, 0, 0, lc} ;
Point(4) = {rr, l, 0, lc} ;
Point(5) = {r, ls, 0, lc} ;
Point(6) = {0, l, 0, lc} ;
Point(7) = {-r,ls , 0, lc} ;
Point(9) = {r,ls+l , 0, lc} ;
Point(10) = {-r,ls+l , 0, lc} ;
Point(11) = {0,ls+ls , 0, lc} ;
/////

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {6, 7};
Line(8) = {7, 10};
Line(9) = {10, 11};
Line(10) = {11, 9};
Line(11) = {9, 5};


Line Loop(2) = {9, 10, 11, 5, 7, 8};
Plane Surface(2) = {2};
Line Loop(1) = {-3, -4,- 5,- 6, -1,- 2};
Plane Surface(1) = {1};

//

For t In {1:20}
  rr += 14.7  ;
my_new_surfs[] = Translate {rr-14.7 , 0, 0} { Duplicata{ Surface{1}; }  };
EndFor


For t In {1:21} //  13/11/9 column
  rr1 += 14.7  ;
my_new_surfs[] = Translate {rr1-29.4,ll, 0} { Duplicata{ Surface{1}; }  };
EndFor

For t In {1:21}
  rr2 += 14.7  ;
my_new_surfs[] = Translate {rr2-29.4,-ll, 0} { Duplicata{ Surface{1}; }  };
EndFor

For t In {1:21}
  rr3 += 14.7  ;
my_new_surfs[] = Translate {rr3-29.4 ,ll*2, 0} { Duplicata{ Surface{1}; }  };
EndFor

//.......
For t In {1:21}
  rr4 += 14.7  ;
my_new_surfs[] = Translate {rr4-29.4 ,-ll*2, 0} { Duplicata{ Surface{1}; }  };
EndFor

/////
For t In {1:19}
  rr5 += 14.7  ;
my_new_surfs[] = Translate {rr5-14.7,ll*3, 0} { Duplicata{ Surface{1}; }  };
EndFor
//.......
For t In {1:19}
  rr6 += 14.7  ;
my_new_surfs[] = Translate {rr6-14.7,-ll*3, 0} { Duplicata{ Surface{1}; }  };
EndFor

//.......
For t In {1:17}
  rr7 += 14.7  ;
my_new_surfs[] = Translate {rr7-0,ll*4, 0} { Duplicata{ Surface{1}; }  };
EndFor
//.......
For t In {1:17}
  rr8 += 14.7  ;
my_new_surfs[] = Translate {rr8-0,-ll*4, 0} { Duplicata{ Surface{1}; }  };
EndFor

//
For t In {1:13}
  rr9 += 14.7  ;
my_new_surfs[] = Translate {rr9+29.4,ll*5, 0} { Duplicata{ Surface{1}; }  };
EndFor
//.......
For t In {1:13}
  rr10 += 14.7  ;
my_new_surfs[] = Translate {rr10+29.4,-ll*5, 0} { Duplicata{ Surface{1}; }  };
EndFor

For t In {1:5}
  rr11 += 14.7  ;
my_new_surfs[] = Translate {rr11+88.2,ll*6, 0} { Duplicata{ Surface{1}; }  };
EndFor
//.......
For t In {1:5}
  rr12 += 14.7  ;
my_new_surfs[] = Translate {rr12+88.2,-ll*6, 0} { Duplicata{ Surface{1}; }  };
EndFor

//.................................................................................................................//



/// .................
For t In {1:21}
  rm += 14.7  ;
my_new_surfs[] = Translate {rm-14.7 , 0, 0} { Duplicata{ Surface{2}; } };

EndFor

For t In {1:22}     //14/12/10/6
  rmm += 14.7  ;
my_new_surfs[] = Translate {rmm -29.4,-ll, 0} { Duplicata{ Surface{2}; }  };
EndFor

For t In {1:20}
  rm2 += 14.7  ;
my_new_surfs[] = Translate {rm2 -14.7 ,ll, 0} { Duplicata{ Surface{2}; }  };
EndFor

//.......
For t In {1:20}
  rm3 += 14.7  ;
my_new_surfs[] = Translate {rm3 -14.7 ,-ll*2, 0} { Duplicata{ Surface{2}; }  };
EndFor

//.......
For t In {1:20}
  rm4 += 14.7  ;
my_new_surfs[] = Translate {rm4 -14.7,ll*2, 0} { Duplicata{ Surface{2}; }  };
EndFor
//.......
For t In {1:20}
  rm5 += 14.7  ;
my_new_surfs[] = Translate {rm5 -14.7,-ll*3, 0} { Duplicata{ Surface{2}; }  };
EndFor

//////////
For t In {1:18}
  rm6 += 14.7  ;
my_new_surfs[] = Translate {rm6 -0,ll*3, 0} { Duplicata{ Surface{2}; }  };
EndFor
//.......
For t In {1:18}
  rm7 += 14.7  ;
my_new_surfs[] = Translate {rm7 -0 ,-ll*4, 0} { Duplicata{ Surface{2}; }  };
EndFor

//
For t In {1:14}
  rm8 += 14.7  ;
my_new_surfs[] = Translate {rm8 +29.4,ll*4, 0} { Duplicata{ Surface{2}; }  };
EndFor
//.......
For t In {1:14}
  rm9 += 14.7  ;
my_new_surfs[] = Translate {rm9 +29.4,-ll*5, 0} { Duplicata{ Surface{2}; }  };
EndFor
//
For t In {1:10}
  rm10 += 14.7  ;
my_new_surfs[] = Translate {rm10 +58.8,ll*5, 0} { Duplicata{ Surface{2}; }  };
EndFor
//.......
For t In {1:10}
  rm11 += 14.7  ;
my_new_surfs[] = Translate {rm11 +58.8,-ll*6, 0} { Duplicata{ Surface{2}; }  };
EndFor

Physical Surface(1) = {1426};
Physical Surface(2) = {1433};
Physical Surface(3) = {1440};
Physical Surface(4) = {1447};
Physical Surface(5) = {1454};
Physical Surface(6) = {2276};
Physical Surface(7) = {2283};
Physical Surface(8) = {2287};
Physical Surface(9) = {2291};
Physical Surface(10) = {2295};
Physical Surface(11) = {2299};
Physical Surface(12) = {2303};
Physical Surface(13) = {2307};
Physical Surface(14) = {2311};
Physical Surface(15) = {2315};
Physical Surface(16) = {1244};
Physical Surface(17) = {1251};
Physical Surface(18) = {1258};
Physical Surface(19) = {1265};
Physical Surface(20) = {1272};
Physical Surface(21) = {1279};
Physical Surface(22) = {1286};
Physical Surface(23) = {1293};
Physical Surface(24) = {1300};
Physical Surface(25) = {1307};
Physical Surface(26) = {1314};
Physical Surface(27) = {1321};
Physical Surface(28) = {1328};
//

Physical Surface(29) = {2157};
Physical Surface(30) = {2164};
Physical Surface(31) = {2168};
Physical Surface(32) = {2172};
Physical Surface(33) = {2176};
Physical Surface(34) = {2180};
Physical Surface(35) = {2184};
Physical Surface(36) = {2188};
Physical Surface(37) = {2192};
Physical Surface(38) = {2196};
Physical Surface(39) = {2200};
Physical Surface(40) = {2204};
Physical Surface(41) = {2208};
Physical Surface(42) = {2212};
Physical Surface(43) = {1006};
Physical Surface(44) = {1013};
Physical Surface(45) = {1020};
Physical Surface(46) = {1027};
Physical Surface(47) = {1034};
Physical Surface(48) = {1041};
Physical Surface(49) = {1048};
Physical Surface(50) = {1055};
Physical Surface(51) = {1062};
Physical Surface(52) = {1069};
Physical Surface(53) = {1076};
Physical Surface(54) = {1083};
Physical Surface(55) = {1090};
Physical Surface(56) = {1097};
Physical Surface(57) = {1104};
Physical Surface(58) = {1111};
Physical Surface(59) = {1118};
Physical Surface(60) = {2006};
Physical Surface(61) = {2013};
Physical Surface(62) = {2017};
Physical Surface(63) = {2021};
Physical Surface(64) = {2025};
Physical Surface(65) = {2029};
Physical Surface(66) = {2033};
Physical Surface(67) = {2037};
Physical Surface(68) = {2041};
Physical Surface(69) = {2045};
Physical Surface(70) = {2049};
Physical Surface(71) = {2053};
Physical Surface(72) = {2057};
Physical Surface(73) = {2061};
Physical Surface(74) = {2065};
Physical Surface(75) = {2069};
Physical Surface(76) = {2073};
Physical Surface(77) = {2077};
Physical Surface(78) = {740};
Physical Surface(79) = {747};
Physical Surface(80) = {754};
Physical Surface(81) = {761};
Physical Surface(82) = {768};
Physical Surface(83) = {775};
Physical Surface(84) = {782};
Physical Surface(85) = {789};
Physical Surface(86) = {796};
Physical Surface(87) = {803};
Physical Surface(88) = {810};
Physical Surface(89) = {817};
Physical Surface(90) = {824};
Physical Surface(91) = {831};
Physical Surface(92) = {838};
Physical Surface(93) = {845};
Physical Surface(94) = {852};
Physical Surface(95) = {859};
Physical Surface(96) = {866};
Physical Surface(97) = {1839};
Physical Surface(98) = {1846};
Physical Surface(99) = {1850};
Physical Surface(100) = {1854};
Physical Surface(101) = {1858};
Physical Surface(102) = {1862};
Physical Surface(103) = {1866};
Physical Surface(104) = {1870};
Physical Surface(105) = {1874};
Physical Surface(106) = {1878};
Physical Surface(107) = {1882};
Physical Surface(108) = {1886};
Physical Surface(109) = {1890};
Physical Surface(110) = {1894};
Physical Surface(111) = {1898};
Physical Surface(112) = {1902};
Physical Surface(113) = {1906};
Physical Surface(114) = {1910};
Physical Surface(115) = {1914};
Physical Surface(116) = {1918};
Physical Surface(117) = {446};
Physical Surface(118) = {453};
Physical Surface(119) = {460};
Physical Surface(120) = {467};
Physical Surface(121) = {474};
Physical Surface(122) = {481};
Physical Surface(123) = {488};
Physical Surface(124) = {495};
Physical Surface(125) = {502};
Physical Surface(126) = {509};
Physical Surface(127) = {516};
Physical Surface(128) = {523};
Physical Surface(129) = {530};
Physical Surface(130) = {537};
Physical Surface(131) = {544};
Physical Surface(132) = {551};
Physical Surface(133) = {558};
Physical Surface(134) = {565};
Physical Surface(135) = {572};
Physical Surface(136) = {579};
Physical Surface(137) = {586};
Physical Surface(138) = {1673};
Physical Surface(139) = {1680};
Physical Surface(140) = {1684};
Physical Surface(141) = {1688};
Physical Surface(142) = {1692};
Physical Surface(143) = {1696};
Physical Surface(144) = {1700};
Physical Surface(145) = {1704};
Physical Surface(146) = {1708};
Physical Surface(147) = {1712};
Physical Surface(148) = {1716};
Physical Surface(149) = {1720};
Physical Surface(150) = {1724};
Physical Surface(151) = {1728};
Physical Surface(152) = {1732};
Physical Surface(153) = {1736};
Physical Surface(154) = {1740};
Physical Surface(155) = {1744};
Physical Surface(156) = {1748};
Physical Surface(157) = {1752};
Physical Surface(158) = {152};
Physical Surface(159) = {159};
Physical Surface(160) = {166};
Physical Surface(161) = {173};
Physical Surface(162) = {180};
Physical Surface(163) = {187};
Physical Surface(164) = {194};
Physical Surface(165) = {201};
Physical Surface(166) = {208};
Physical Surface(167) = {215};
Physical Surface(168) = {222};
Physical Surface(169) = {229};
Physical Surface(170) = {236};
Physical Surface(171) = {243};
Physical Surface(172) = {250};
Physical Surface(173) = {257};
Physical Surface(174) = {264};
Physical Surface(175) = {271};
Physical Surface(176) = {278};
Physical Surface(177) = {285};
Physical Surface(178) = {292};
Physical Surface(179) = {2};
Physical Surface(180) = {1496};
Physical Surface(181) = {1500};
Physical Surface(182) = {1504};
Physical Surface(183) = {1508};
Physical Surface(184) = {1512};
Physical Surface(185) = {1516};
Physical Surface(186) = {1520};
Physical Surface(187) = {1524};
Physical Surface(188) = {1528};
Physical Surface(189) = {1532};
Physical Surface(190) = {1536};
Physical Surface(191) = {1540};
Physical Surface(192) = {1544};
Physical Surface(193) = {1548};
Physical Surface(194) = {1552};
Physical Surface(195) = {1556};
Physical Surface(196) = {1560};
Physical Surface(197) = {1564};
Physical Surface(198) = {1568};
Physical Surface(199) = {1572};
Physical Surface(200) = {1576};
Physical Surface(201) = {1};
Physical Surface(202) = {12};
Physical Surface(203) = {19};
Physical Surface(204) = {26};
Physical Surface(205) = {33};
Physical Surface(206) = {40};
Physical Surface(207) = {47};
Physical Surface(208) = {54};
Physical Surface(209) = {61};
Physical Surface(210) = {68};
Physical Surface(211) = {75};
Physical Surface(212) = {82};
Physical Surface(213) = {89};
Physical Surface(214) = {96};
Physical Surface(215) = {103};
Physical Surface(216) = {110};
Physical Surface(217) = {117};
Physical Surface(218) = {124};
Physical Surface(219) = {131};
Physical Surface(220) = {138};
Physical Surface(221) = {145};
Physical Surface(222) = {1581};
Physical Surface(223) = {1588};
Physical Surface(224) = {1592};
Physical Surface(225) = {1596};
Physical Surface(226) = {1600};
Physical Surface(227) = {1604};
Physical Surface(228) = {1608};
Physical Surface(229) = {1612};
Physical Surface(230) = {1616};
Physical Surface(231) = {1620};
Physical Surface(232) = {1624};
Physical Surface(233) = {1628};
Physical Surface(234) = {1632};
Physical Surface(235) = {1636};
Physical Surface(236) = {1640};
Physical Surface(237) = {1644};
Physical Surface(238) = {1648};
Physical Surface(239) = {1652};
Physical Surface(240) = {1656};
Physical Surface(241) = {1660};
Physical Surface(242) = {1664};
Physical Surface(243) = {1668};
Physical Surface(244) = {299};
Physical Surface(245) = {306};
Physical Surface(246) = {313};
Physical Surface(247) = {320};
Physical Surface(248) = {327};
Physical Surface(249) = {334};
Physical Surface(250) = {341};
Physical Surface(251) = {348};
Physical Surface(252) = {355};
Physical Surface(253) = {362};
Physical Surface(254) = {369};
Physical Surface(255) = {376};
Physical Surface(256) = {383};
Physical Surface(257) = {390};
Physical Surface(258) = {397};
Physical Surface(259) = {404};
Physical Surface(260) = {411};
Physical Surface(261) = {418};
Physical Surface(262) = {425};
Physical Surface(263) = {432};
Physical Surface(264) = {439};
Physical Surface(265) = {1756};
Physical Surface(266) = {1763};
Physical Surface(267) = {1767};
Physical Surface(268) = {1771};
Physical Surface(269) = {1775};
Physical Surface(270) = {1779};
Physical Surface(271) = {1783};
Physical Surface(272) = {1787};
Physical Surface(273) = {1791};
Physical Surface(274) = {1795};
Physical Surface(275) = {1799};
Physical Surface(276) = {1803};
Physical Surface(277) = {1807};
Physical Surface(278) = {1811};
Physical Surface(279) = {1815};
Physical Surface(280) = {1819};
Physical Surface(281) = {1823};
Physical Surface(282) = {1827};
Physical Surface(283) = {1831};
Physical Surface(284) = {1835};
Physical Surface(285) = {593};
Physical Surface(286) = {600};
Physical Surface(287) = {607};
Physical Surface(288) = {614};
Physical Surface(289) = {621};
Physical Surface(290) = {628};
Physical Surface(291) = {635};
Physical Surface(292) = {642};
Physical Surface(293) = {649};
Physical Surface(294) = {656};
Physical Surface(295) = {663};
Physical Surface(296) = {670};
Physical Surface(297) = {677};
Physical Surface(298) = {684};
Physical Surface(299) = {691};
Physical Surface(300) = {698};
Physical Surface(301) = {705};
Physical Surface(302) = {712};
Physical Surface(303) = {719};
Physical Surface(304) = {726};
Physical Surface(305) = {733};
Physical Surface(306) = {1922};
Physical Surface(307) = {1929};
Physical Surface(308) = {1933};
Physical Surface(309) = {1937};
Physical Surface(310) = {1941};
Physical Surface(311) = {1945};
Physical Surface(312) = {1949};
Physical Surface(313) = {1953};
Physical Surface(314) = {1957};
Physical Surface(315) = {1961};
Physical Surface(316) = {1965};
Physical Surface(317) = {1969};
Physical Surface(318) = {1973};
Physical Surface(319) = {1977};
Physical Surface(320) = {1981};
Physical Surface(321) = {1985};
Physical Surface(322) = {1989};
Physical Surface(323) = {1993};
Physical Surface(324) = {1997};
Physical Surface(325) = {2001};
Physical Surface(326) = {873};
Physical Surface(327) = {880};
Physical Surface(328) = {887};
Physical Surface(329) = {894};
Physical Surface(330) = {901};
Physical Surface(331) = {908};
Physical Surface(332) = {915};
Physical Surface(333) = {922};
Physical Surface(334) = {929};
Physical Surface(335) = {936};
Physical Surface(336) = {943};
Physical Surface(337) = {950};
Physical Surface(338) = {957};
Physical Surface(339) = {964};
Physical Surface(340) = {971};
Physical Surface(341) = {978};
Physical Surface(342) = {985};
Physical Surface(343) = {992};
Physical Surface(344) = {999};
Physical Surface(345) = {2081};
Physical Surface(346) = {2088};
Physical Surface(347) = {2092};
Physical Surface(348) = {2096};
Physical Surface(349) = {2100};
Physical Surface(350) = {2104};
Physical Surface(351) = {2108};
Physical Surface(352) = {2112};
Physical Surface(353) = {2116};
Physical Surface(354) = {2120};
Physical Surface(355) = {2124};
Physical Surface(356) = {2128};
Physical Surface(357) = {2132};
Physical Surface(358) = {2136};
Physical Surface(359) = {2140};
Physical Surface(360) = {2144};
Physical Surface(361) = {2148};
Physical Surface(362) = {2152};
Physical Surface(363) = {1125};
Physical Surface(364) = {1132};
Physical Surface(365) = {1139};
Physical Surface(366) = {1146};
Physical Surface(367) = {1153};
Physical Surface(368) = {1160};
Physical Surface(369) = {1167};
Physical Surface(370) = {1174};
Physical Surface(371) = {1181};
Physical Surface(372) = {1188};
Physical Surface(373) = {1195};
Physical Surface(374) = {1202};
Physical Surface(375) = {1209};
Physical Surface(376) = {1216};
Physical Surface(377) = {1223};
Physical Surface(378) = {1230};
Physical Surface(379) = {1237};
Physical Surface(380) = {2216};
Physical Surface(381) = {2223};
Physical Surface(382) = {2227};
Physical Surface(383) = {2231};
Physical Surface(384) = {2235};
Physical Surface(385) = {2239};
Physical Surface(386) = {2243};
Physical Surface(387) = {2247};
Physical Surface(388) = {2251};
Physical Surface(389) = {2255};
Physical Surface(390) = {2259};
Physical Surface(391) = {2263};
Physical Surface(392) = {2267};
Physical Surface(393) = {2271};
Physical Surface(394) = {1335};
Physical Surface(395) = {1342};
Physical Surface(396) = {1349};
Physical Surface(397) = {1356};
Physical Surface(398) = {1363};
Physical Surface(399) = {1370};
Physical Surface(400) = {1377};
Physical Surface(401) = {1384};
Physical Surface(402) = {1391};
Physical Surface(403) = {1398};
Physical Surface(404) = {1405};
Physical Surface(405) = {1412};
Physical Surface(406) = {1419};
Physical Surface(407) = {2319};
Physical Surface(408) = {2326};
Physical Surface(409) = {2332};
Physical Surface(410) = {2338};
Physical Surface(411) = {2342};
Physical Surface(412) = {2346};
Physical Surface(413) = {2350};
Physical Surface(414) = {2354};
Physical Surface(415) = {2359};
Physical Surface(416) = {2365};
Physical Surface(417) = {1461};
Physical Surface(418) = {1468};
Physical Surface(419) = {1475};
Physical Surface(420) = {1482};
Physical Surface(421) = {1489};
Physical Line(1) = {6, 1582, 1587, 1586, 303, 302, 1762, 598, 597, 596, 1928, 1927, 877, 876, 2087, 2086, 1129, 1128, 1127, 1135, 2222, 2221, 1339, 1338, 1337, 1345, 2325, 2324, 2323, 2331, 2330, 2337, 1465, 1464, 1463, 1471, 1470, 1478, 1477, 1485, 1484, 1492, 1491, 1490, 2358, 2364, 2363, 2370, 2369, 2368, 1414, 1422, 1421, 1420, 2275, 2274, 1232, 1240, 1239, 1238, 2156, 2155, 1001, 1000, 2005, 2004, 735, 734, 739, 1838, 441, 440, 1672, 1671, 1670, 146, 1580, 1579, 1578, 293, 298, 1755, 588, 587, 592, 1921, 1920, 867, 872, 2080, 2079, 1119, 1124, 1123, 1117, 2215, 2214, 1329, 1334, 1333, 1327, 2318, 2317, 2316, 2313, 2312, 2309, 1455, 1460, 1459, 1453, 1452, 1446, 1445, 1439, 1438, 1432, 1431, 1430, 2288, 2285, 2284, 2278, 2277, 2282, 1256, 1250, 1249, 1248, 2158, 2163, 1018, 1012, 1011, 1010, 2007, 2012, 745, 744, 1840, 1845, 451, 450, 449, 1679, 157, 156, 9, 8, 7};
