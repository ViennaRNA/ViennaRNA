#ifndef CONFIG_ARC_STRUCT_H
#define CONFIG_ARC_STRUCT_H

/**
 * @brief The configArc struct holds the description of a single arc of the config.
 */
typedef struct
{
    /// description of the arc
    /// - number of arc segments
    int numberOfArcSegments;

    /// - arcAngle
    double arcAngle;
} configArc;

#endif
