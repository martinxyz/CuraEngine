/** Copyright (C) 2013 David Braam - Released under terms of the AGPLv3 License */
#include <stdarg.h>

#include "gcodeExport.h"
#include "pathOrderOptimizer.h"
#include "timeEstimate.h"
#include "settings.h"
#include "utils/logoutput.h"

#if defined(__APPLE__) && defined(__MACH__)
//On MacOS the file offset functions are always 64bit.
#define off64_t off_t
#define ftello64 ftello
#define fseeko64 fseeko
#endif

GCodeExport::GCodeExport()
: currentPosition(0,0,0)
{
    extrusionAmount = 0;
    extrusionPerMM = 0;
    retractionAmount = 4.5;
    minimalExtrusionBeforeRetraction = 0.0;
    extrusionAmountAtPreviousRetraction = -10000;
    extruderSwitchRetraction = 14.5;
    extruderNr = 0;
    currentFanSpeed = -1;
    
    totalPrintTime = 0.0;
    for(unsigned int e=0; e<MAX_EXTRUDERS; e++)
        totalFilament[e] = 0.0;
    
    currentSpeed = 0.0;
    retractionSpeed = 45;
    isRetracted = true;
    setFlavor(GCODE_FLAVOR_REPRAP);
    memset(extruderOffset, 0, sizeof(extruderOffset));
    f = stdout;
}

GCodeExport::~GCodeExport()
{
    if (f && f != stdout)
        fclose(f);
}

void GCodeExport::replaceTagInStart(const char* tag, const char* replaceValue)
{
    if (f == stdout)
    {
        log("Replace:%s:%s\n", tag, replaceValue);
        return;
    }
    off64_t oldPos = ftello64(f);
    
    char buffer[1024];
    fseeko64(f, 0, SEEK_SET);
    fread(buffer, 1024, 1, f);
    
    char* c = strstr(buffer, tag);
    memset(c, ' ', strlen(tag));
    if (c) memcpy(c, replaceValue, strlen(replaceValue));
    
    fseeko64(f, 0, SEEK_SET);
    fwrite(buffer, 1024, 1, f);
    
    fseeko64(f, oldPos, SEEK_SET);
}

void GCodeExport::setExtruderOffset(int id, Point p)
{
    extruderOffset[id] = p;
}

void GCodeExport::setFlavor(int flavor)
{
    this->flavor = flavor;
    if (flavor == GCODE_FLAVOR_MACH3)
        for(int n=0; n<MAX_EXTRUDERS; n++)
            extruderCharacter[n] = 'A' + n;
    else
        for(int n=0; n<MAX_EXTRUDERS; n++)
            extruderCharacter[n] = 'E';
}
int GCodeExport::getFlavor()
{
    return this->flavor;
}

void GCodeExport::setFilename(const char* filename)
{
    f = fopen(filename, "w+");
}

bool GCodeExport::isOpened()
{
    return f != NULL;
}

void GCodeExport::setExtrusion(int layerThickness, int filamentDiameter, int flow)
{
    double filamentArea = M_PI * (INT2MM(filamentDiameter) / 2.0) * (INT2MM(filamentDiameter) / 2.0);
    if (flavor == GCODE_FLAVOR_ULTIGCODE)//UltiGCode uses volume extrusion as E value, and thus does not need the filamentArea in the mix.
        extrusionPerMM = INT2MM(layerThickness);
    else
        extrusionPerMM = INT2MM(layerThickness) / filamentArea * double(flow) / 100.0;
}

void GCodeExport::setRetractionSettings(int retractionAmount, int retractionSpeed, int extruderSwitchRetraction, int minimalExtrusionBeforeRetraction, int zHop, int retractionAmountPrime)
{
    this->retractionAmount = INT2MM(retractionAmount);
    this->retractionAmountPrime = INT2MM(retractionAmountPrime);
    this->retractionSpeed = retractionSpeed;
    this->extruderSwitchRetraction = INT2MM(extruderSwitchRetraction);
    this->minimalExtrusionBeforeRetraction = INT2MM(minimalExtrusionBeforeRetraction);
    this->retractionZHop = zHop;
}

void GCodeExport::setZ(int z)
{
    this->zPos = z;
}

Point GCodeExport::getPositionXY()
{
    return Point(currentPosition.x, currentPosition.y);
}

int GCodeExport::getPositionZ()
{
    return currentPosition.z;
}

int GCodeExport::getExtruderNr()
{
    return extruderNr;
}

double GCodeExport::getTotalFilamentUsed(int e)
{
    if (e == extruderNr)
        return totalFilament[e] + extrusionAmount;
    return totalFilament[e];
}

double GCodeExport::getTotalPrintTime()
{
    return totalPrintTime;
}

void GCodeExport::updateTotalPrintTime()
{
    totalPrintTime += estimateCalculator.calculate();
    estimateCalculator.reset();
}

void GCodeExport::writeComment(const char* comment, ...)
{
    va_list args;
    va_start(args, comment);
    fprintf(f, ";");
    vfprintf(f, comment, args);
    fprintf(f, "\n");
    va_end(args);
}

void GCodeExport::writeLine(const char* line, ...)
{
    va_list args;
    va_start(args, line);
    vfprintf(f, line, args);
    fprintf(f, "\n");
    va_end(args);
}

void GCodeExport::resetExtrusionValue()
{
    if (extrusionAmount != 0.0 && flavor != GCODE_FLAVOR_MAKERBOT && flavor != GCODE_FLAVOR_BFB)
    {
        fprintf(f, "G92 %c0\n", extruderCharacter[extruderNr]);
        totalFilament[extruderNr] += extrusionAmount;
        extrusionAmountAtPreviousRetraction -= extrusionAmount;
        extrusionAmount = 0.0;
    }
}

void GCodeExport::writeDelay(double timeAmount)
{
    fprintf(f, "G4 P%d\n", int(timeAmount * 1000));
    totalPrintTime += timeAmount;
}

void GCodeExport::writeMove(Point p, double speed, int lineWidth)
{
    if (flavor == GCODE_FLAVOR_BFB)
    {
        //For Bits From Bytes machines, we need to handle this completely differently. As they do not use E values but RPM values.
        double fspeed = speed * 60;
        double rpm = (extrusionPerMM * double(lineWidth) / 1000.0) * speed * 60;
        const double mm_per_rpm = 4.0; //All BFB machines have 4mm per RPM extrusion.
        rpm /= mm_per_rpm;
        if (rpm > 0)
        {
            if (isRetracted)
            {
                if (currentSpeed != int(rpm * 10))
                {
                    //fprintf(f, "; %f e-per-mm %d mm-width %d mm/s\n", extrusionPerMM, lineWidth, speed);
                    fprintf(f, "M108 S%0.1f\n", rpm);
                    currentSpeed = int(rpm * 10);
                }
                fprintf(f, "M%d01\n", extruderNr);
                isRetracted = false;
            }
            //Fix the speed by the actual RPM we are asking, because of rounding errors we cannot get all RPM values, but we have a lot more resolution in the feedrate value.
            // (Trick copied from KISSlicer, thanks Jonathan)
            fspeed *= (rpm / (roundf(rpm * 100) / 100));

            //Increase the extrusion amount to calculate the amount of filament used.
            Point diff = p - getPositionXY();
            extrusionAmount += extrusionPerMM * INT2MM(lineWidth) * vSizeMM(diff);
        }else{
            //If we are not extruding, check if we still need to disable the extruder. This causes a retraction due to auto-retraction.
            if (!isRetracted)
            {
                fprintf(f, "M103\n");
                isRetracted = true;
            }
        }
        fprintf(f, "G1 X%0.2f Y%0.2f Z%0.2f F%0.1f\n", INT2MM(p.X - extruderOffset[extruderNr].X), INT2MM(p.Y - extruderOffset[extruderNr].Y), INT2MM(zPos), fspeed);
    }else{
        
        //Normal E handling.
        if (lineWidth != 0)
        {
            Point diff = p - getPositionXY();
            if (isRetracted)
            {
                if (retractionZHop > 0)
                    fprintf(f, "G1 Z%0.2f\n", double(currentPosition.z)/1000);
                if (flavor == GCODE_FLAVOR_ULTIGCODE)
                {
                    fprintf(f, "G11\n");
                }else{
                    extrusionAmount += retractionAmountPrime;
                    fprintf(f, "G1 F%i %c%0.5lf\n", retractionSpeed * 60, extruderCharacter[extruderNr], extrusionAmount);
                    currentSpeed = retractionSpeed;
                    estimateCalculator.plan(TimeEstimateCalculator::Position(INT2MM(p.X), INT2MM(p.Y), INT2MM(zPos), extrusionAmount), currentSpeed);
                }
                if (extrusionAmount > 10000.0) //According to https://github.com/Ultimaker/CuraEngine/issues/14 having more then 21m of extrusion causes inaccuracies. So reset it every 10m, just to be sure.
                    resetExtrusionValue();
                isRetracted = false;
            }
            extrusionAmount += extrusionPerMM * INT2MM(lineWidth) * vSizeMM(diff);
            fprintf(f, "G1");
        }else{
            fprintf(f, "G0");
        }

        if (currentSpeed != speed)
        {
            fprintf(f, " F%i", int(speed * 60));
            currentSpeed = speed;
        }

        fprintf(f, " X%0.2f Y%0.2f", INT2MM(p.X - extruderOffset[extruderNr].X), INT2MM(p.Y - extruderOffset[extruderNr].Y));
        if (zPos != currentPosition.z)
            fprintf(f, " Z%0.2f", INT2MM(zPos));
        if (lineWidth != 0)
            fprintf(f, " %c%0.5lf", extruderCharacter[extruderNr], extrusionAmount);
        fprintf(f, "\n");
    }
    
    currentPosition = Point3(p.X, p.Y, zPos);
    estimateCalculator.plan(TimeEstimateCalculator::Position(INT2MM(currentPosition.x), INT2MM(currentPosition.y), INT2MM(currentPosition.z), extrusionAmount), speed);
}

void GCodeExport::writeRetraction()
{
    if (flavor == GCODE_FLAVOR_BFB)//BitsFromBytes does automatic retraction.
        return;
    
    if (retractionAmount > 0 && !isRetracted && extrusionAmountAtPreviousRetraction + minimalExtrusionBeforeRetraction < extrusionAmount)
    {
        if (flavor == GCODE_FLAVOR_ULTIGCODE)
        {
            fprintf(f, "G10\n");
        }else{
            fprintf(f, "G1 F%i %c%0.5lf\n", retractionSpeed * 60, extruderCharacter[extruderNr], extrusionAmount - retractionAmount);
            currentSpeed = retractionSpeed;
            estimateCalculator.plan(TimeEstimateCalculator::Position(INT2MM(currentPosition.x), INT2MM(currentPosition.y), INT2MM(currentPosition.z), extrusionAmount - retractionAmount), currentSpeed);
        }
        if (retractionZHop > 0)
            fprintf(f, "G1 Z%0.2f\n", INT2MM(currentPosition.z + retractionZHop));
        extrusionAmountAtPreviousRetraction = extrusionAmount;
        isRetracted = true;
    }
}

void GCodeExport::switchExtruder(int newExtruder)
{
    if (extruderNr == newExtruder)
        return;
    if (flavor == GCODE_FLAVOR_BFB)
    {
        if (!isRetracted)
            fprintf(f, "M103\n");
        isRetracted = true;
        return;
    }
    
    if (flavor == GCODE_FLAVOR_ULTIGCODE)
    {
        fprintf(f, "G10 S1\n");
    }else{
        fprintf(f, "G1 F%i %c%0.5lf\n", retractionSpeed * 60, extruderCharacter[extruderNr], extrusionAmount - extruderSwitchRetraction);
        currentSpeed = retractionSpeed;
    }
    resetExtrusionValue();
    extruderNr = newExtruder;
    if (flavor == GCODE_FLAVOR_MACH3)
        resetExtrusionValue();
    isRetracted = true;
    if (flavor == GCODE_FLAVOR_MAKERBOT)
        fprintf(f, "M135 T%i\n", extruderNr);
    else
        fprintf(f, "T%i\n", extruderNr);
}

void GCodeExport::writeCode(const char* str)
{
    fprintf(f, "%s\n", str);
}

void GCodeExport::writeFanCommand(int speed)
{
    if (currentFanSpeed == speed)
        return;
    if (speed > 0)
    {
        if (flavor == GCODE_FLAVOR_MAKERBOT)
            fprintf(f, "M126 T0 ; value = %d\n", speed * 255 / 100);
        else
            fprintf(f, "M106 S%d\n", speed * 255 / 100);
    }
    else
    {
        if (flavor == GCODE_FLAVOR_MAKERBOT)
            fprintf(f, "M127 T0\n");
        else
            fprintf(f, "M107\n");
    }
    currentFanSpeed = speed;
}

int GCodeExport::getFileSize(){
    return ftell(f);
}
void GCodeExport::tellFileSize() {
    float fsize = (float) ftell(f);
    if(fsize > 1024*1024) {
        fsize /= 1024.0*1024.0;
        log("Wrote %5.1f MB.\n",fsize);
    }
    if(fsize > 1024) {
        fsize /= 1024.0;
        log("Wrote %5.1f kilobytes.\n",fsize);
    }
}

void GCodeExport::finalize(int maxObjectHeight, int moveSpeed, const char* endCode)
{
    writeFanCommand(0);
    writeRetraction();
    setZ(maxObjectHeight + 5000);
    writeMove(getPositionXY(), moveSpeed, 0);
    writeCode(endCode);
    log("Print time: %d\n", int(getTotalPrintTime()));
    log("Filament: %d\n", int(getTotalFilamentUsed(0)));
    log("Filament2: %d\n", int(getTotalFilamentUsed(1)));
    
    if (getFlavor() == GCODE_FLAVOR_ULTIGCODE)
    {
        char numberString[16];
        sprintf(numberString, "%d", int(getTotalPrintTime()));
        replaceTagInStart("<__TIME__>", numberString);
        sprintf(numberString, "%d", int(getTotalFilamentUsed(0)));
        replaceTagInStart("<FILAMENT>", numberString);
        sprintf(numberString, "%d", int(getTotalFilamentUsed(1)));
        replaceTagInStart("<FILAMEN2>", numberString);
    }
}

GCodePath* GCodePlanner::getLatestPathWithConfig(GCodePathConfig* config)
{
    if (paths.size() > 0 && paths[paths.size()-1].config == config && !paths[paths.size()-1].done)
        return &paths[paths.size()-1];
    paths.push_back(GCodePath());
    GCodePath* ret = &paths[paths.size()-1];
    ret->retract = false;
    ret->config = config;
    ret->extruder = currentExtruder;
    ret->done = false;
    return ret;
}

void GCodePlanner::addPointWithConfig(Point p, GCodePathConfig* config)
{
    GCodePath *path = getLatestPathWithConfig(config);
    path->segments.push_back(GCodePathSegment(lastPosition, p, config));
    lastPosition = p;
}

void GCodePlanner::forceNewPathStart()
{
    if (paths.size() > 0)
        paths[paths.size()-1].done = true;
}

GCodePlanner::GCodePlanner(Point startPositionXY, int startExtruder, int travelSpeed, int retractionMinimalDistance)
: travelConfig(travelSpeed, 0, "travel")
{
    startPosition = startPositionXY;
    lastPosition = startPositionXY;
    comb = NULL;
    extraTime = 0.0;
    layerTime = 0.0;
    layerTimeWithoutSpeedLimits = 0.0;
    forceRetraction = false;
    alwaysRetract = false;
    currentExtruder = startExtruder;
    this->retractionMinimalDistance = retractionMinimalDistance;
}
GCodePlanner::~GCodePlanner()
{
    if (comb)
        delete comb;
}

void GCodePlanner::addTravel(Point p)
{
    GCodePath* path = getLatestPathWithConfig(&travelConfig);
    if (forceRetraction)
    {
        if (!shorterThen(lastPosition - p, retractionMinimalDistance))
        {
            path->retract = true;
        }
        forceRetraction = false;
    }else if (comb != NULL)
    {
        vector<Point> pointList;
        if (comb->calc(lastPosition, p, pointList))
        {
            for(unsigned int n=0; n<pointList.size(); n++)
            {
                addPointWithConfig(p, &travelConfig);
            }
        }else{
            if (!shorterThen(lastPosition - p, retractionMinimalDistance))
                path->retract = true;
        }
    }else if (alwaysRetract)
    {
        if (!shorterThen(lastPosition - p, retractionMinimalDistance))
            path->retract = true;
    }
    addPointWithConfig(p, &travelConfig);
}

void GCodePlanner::moveInsideCombBoundary(int distance)
{
    if (!comb || comb->checkInside(lastPosition)) return;
    Point p = lastPosition;
    if (comb->moveInside(&p, distance))
    {
        //Move inside again, so we move out of tight 90deg corners
        comb->moveInside(&p, distance);
        if (comb->checkInside(p))
        {
            addTravel(p);
            //Make sure the that any retraction happens after this move, not before it by starting a new move path.
            forceNewPathStart();
        }
    }
}

void GCodePlanner::addPolygon(PolygonRef polygon, int startIdx, GCodePathConfig* config)
{
    Point p0 = polygon[startIdx];
    addTravel(p0);
    for(unsigned int i=1; i<polygon.size(); i++)
    {
        Point p1 = polygon[(startIdx + i) % polygon.size()];
        addPointWithConfig(p1, config);
        p0 = p1;
    }
    if (polygon.size() > 2)
        addPointWithConfig(polygon[startIdx], config);
}

void GCodePlanner::addPolygonsByOptimizer(Polygons& polygons, GCodePathConfig* config)
{
    PathOrderOptimizer orderOptimizer(lastPosition);
    for(unsigned int i=0;i<polygons.size();i++)
        orderOptimizer.addPolygon(polygons[i]);
    orderOptimizer.optimize();
    for(unsigned int i=0;i<orderOptimizer.polyOrder.size();i++)
    {
        int nr = orderOptimizer.polyOrder[i];
        addPolygon(polygons[nr], orderOptimizer.polyStart[nr], config);
    }
}

void GCodePlanner::simpleTimeEstimate(double &travelTime, double &extrudeTime)
{
    travelTime = 0.0;
    extrudeTime = 0.0;
    for(unsigned int n=0; n<paths.size(); n++)
    {
        GCodePath* path = &paths[n];
        for(unsigned int i=0; i<path->segments.size(); i++)
        {
            GCodePathSegment &seg = path->segments[i];
            double thisTime = seg.length / seg.speed;
            if (seg.lineWidth != 0)
                extrudeTime += thisTime;
            else
                travelTime += thisTime;
        }
    }
}

void GCodePlanner::enforceSpeedLimits(double minTime, int minSpeed, int maxSpeed)
{
    double travelTime, extrudeTime, totalTime;
    simpleTimeEstimate(travelTime, extrudeTime);
    totalTime = extrudeTime + travelTime;
    this->layerTimeWithoutSpeedLimits = totalTime;

    double factor = 1.0;
    if (totalTime < minTime && extrudeTime > 0.0)
    {
        double minExtrudeTime = minTime - travelTime;
        if (minExtrudeTime < 1)
            minExtrudeTime = 1;
        factor = extrudeTime / minExtrudeTime;
    }

    // OPTIMIZE: check if significant time is lost when doing this on huge layers where it is not needed. Also, integer-float conversion could be eliminated by making "factor" integer.
    for(unsigned int n=0; n<paths.size(); n++)
    {
        for(unsigned int i=0; i<paths[n].segments.size(); i++)
        {
            GCodePathSegment &seg = paths[n].segments[i];
            if (seg.lineWidth == 0)
                continue;
            double speed = seg.speed * factor;
            if (minSpeed != 0 && speed < minSpeed)
                speed = minSpeed;
            if (maxSpeed != 0 && speed > maxSpeed)
                speed = maxSpeed;
            seg.speed = speed;
        }
    }
        
    simpleTimeEstimate(travelTime, extrudeTime);
    totalTime = extrudeTime + travelTime;

    if (minTime - totalTime > 0.1)
    {
        //TODO: Use up this extra time (circle around the print?)
        this->extraTime = minTime - totalTime;
    }
    this->layerTime = totalTime;
}

double GCodePlanner::limitFlowGrowthRate(double initialFlow2D, double maxFlowGrowthRate, bool forward)
{
    // We limit both decay and growth of the flow.
    // Decay is limited by limiting growth backwards in time.
    double maxFlow = initialFlow2D;
    // OPTIMIZE: check if significant time is lost when doing this on huge layers where it is not needed.
    for(unsigned int n=0; n<paths.size(); n++)
    {
        int nn = forward ? n : (paths.size() - 1 - n);
        for(unsigned int i=0; i<paths[nn].segments.size(); i++)
        {
            int ii = forward ? i : (paths[nn].segments.size() - 1 - i);

            GCodePathSegment &seg = paths[nn].segments[ii];
            if (seg.lineWidth == 0)
                continue;
            double flow = seg.speed * seg.lineWidth;
            if (maxFlow == 0)
                maxFlow = flow;
            if (flow > maxFlow)
            {
                flow = maxFlow;
                seg.speed = flow / seg.lineWidth;
            }
            // We achieve exponential growth (in time) by doing constant growth (per extruded volume).
            maxFlow = flow + seg.length * seg.lineWidth * maxFlowGrowthRate;
        }
    }
    return maxFlow;
}

void GCodePlanner::writeGCode(GCodeExport& gcode, bool liftHeadIfNeeded, int layerThickness)
{
    GCodePathConfig* lastConfig = NULL;
    int extruder = gcode.getExtruderNr();

    for(unsigned int n=0; n<paths.size(); n++)
    {
        GCodePath* path = &paths[n];
        if (extruder != path->extruder)
        {
            extruder = path->extruder;
            gcode.switchExtruder(extruder);
        }else if (path->retract)
        {
            gcode.writeRetraction();
        }
        if (path->config != &travelConfig && lastConfig != path->config)
        {
            gcode.writeComment("TYPE:%s", path->config->name);
            lastConfig = path->config;
        }
        
        if (path->segments.size() == 1 && path->config != &travelConfig && shorterThen(gcode.getPositionXY() - path->segments[0].pos, path->segments[0].lineWidth * 2))
        {
            //Check for lots of small moves and combine them into one large line
            Point p0 = path->segments[0].pos;
            unsigned int i = n + 1;
            while(i < paths.size() && paths[i].segments.size() == 1 && shorterThen(p0 - paths[i].segments[0].pos, path->segments[0].lineWidth * 2))
            {
                p0 = paths[i].segments[0].pos;
                i ++;
            }
            if (paths[i-1].config == &travelConfig)
                i --;
            if (i > n + 2)
            {
                p0 = gcode.getPositionXY();
                for(unsigned int x=n; x<i-1; x+=2)
                {
                    GCodePathSegment s1 = paths[x].segments[0];
                    GCodePathSegment s2 = paths[x+1].segments[0];
                    int64_t oldLen = vSize(p0 - s1.pos);
                    Point newPoint = (s1.pos + s2.pos) / 2;
                    int newSpeed = s1.speed;
                    int64_t newLen = vSize(gcode.getPositionXY() - newPoint);
                    int newLineWidth = s1.lineWidth * oldLen / newLen;
                    if (newLen > 0)
                        gcode.writeMove(newPoint, newSpeed, newLineWidth);
                    
                    p0 = s2.pos;
                }
                GCodePathSegment seg = paths[i-1].segments[0];
                gcode.writeMove(seg.pos, seg.speed, seg.lineWidth);
                n = i - 1;
                continue;
            }
        }
        
        bool spiralize = path->config->spiralize;
        if (spiralize)
        {
            //Check if we are the last spiralize path in the list, if not, do not spiralize.
            for(unsigned int m=n+1; m<paths.size(); m++)
            {
                if (paths[m].config->spiralize)
                    spiralize = false;
            }
        }
        if (spiralize)
        {
            //If we need to spiralize then raise the head slowly by 1 layer as this path progresses.
            float totalLength = 0.0;
            int z = gcode.getPositionZ();
            Point p0 = gcode.getPositionXY();
            for(unsigned int i=0; i<path->segments.size(); i++)
            {
                Point p1 = path->segments[i].pos;
                totalLength += vSizeMM(p0 - p1);
                p0 = p1;
            }
            
            float length = 0.0;
            p0 = gcode.getPositionXY();
            for(unsigned int i=0; i<path->segments.size(); i++)
            {
                GCodePathSegment seg = path->segments[i];
                Point p1 = seg.pos;
                length += vSizeMM(p0 - p1);
                p0 = p1;
                gcode.setZ(z + layerThickness * length / totalLength);
                gcode.writeMove(seg.pos, seg.speed, seg.lineWidth);
            }
        }else{
            for(unsigned int i=0; i<path->segments.size(); i++)
            {
                GCodePathSegment seg = path->segments[i];
                gcode.writeMove(seg.pos, seg.speed, seg.lineWidth);
            }
        }
    }
    
    gcode.updateTotalPrintTime();
    if (liftHeadIfNeeded && extraTime > 0.0)
    {
        gcode.writeComment("Small layer, adding delay of %f", extraTime);
        gcode.writeRetraction();
        gcode.setZ(gcode.getPositionZ() + MM2INT(3.0));
        gcode.writeMove(gcode.getPositionXY(), travelConfig.speed, 0);
        gcode.writeMove(gcode.getPositionXY() - Point(-MM2INT(20.0), 0), travelConfig.speed, 0);
        gcode.writeDelay(extraTime);
    }
}
