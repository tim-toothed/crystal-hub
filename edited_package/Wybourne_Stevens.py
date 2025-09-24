from edited_package.constants import LambdaConstants, LStheta, theta

def WybourneToStevens(ion, Bdict, LS=False):
    StevDict = {}
    for Anm in Bdict:
        n = int(Anm[1])
        m = int(Anm[2:])
        if LS:
            StevDict['B'+Anm[1:]] = LambdaConstants[n][m]*LStheta(ion,n)*Bdict[Anm]
        else:
            StevDict['B'+Anm[1:]] = LambdaConstants[n][m]*theta(ion,n)*Bdict[Anm]
    return StevDict

def StevensToWybourne(ion, Bdict, LS=False):
    WybDict = {}
    for Anm in Bdict:
        n = int(Anm[1])
        m = int(Anm[2:])
        if LS:
            WybDict['B'+Anm[1:]] = Bdict[Anm]/(LambdaConstants[n][m]*LStheta(ion,n))
        else:
            WybDict['B'+Anm[1:]] = Bdict[Anm]/(LambdaConstants[n][m]*theta(ion,n))
    return WybDict