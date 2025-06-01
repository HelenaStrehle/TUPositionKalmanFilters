
function [accX_corr, accY_corr, accZ_corr, GyroX_corr, GyroY_corr, GyroZ_corr] = preProcessSensorData(imuTable, acc_offset, gyro_offset)
    [accX_corr, accY_corr, accZ_corr] = correctAccelerometer(imuTable, acc_offset);
    [GyroX_corr, GyroY_corr, GyroZ_corr] = correctGyroscope(imuTable, gyro_offset);
end

function [accX_corr, accY_corr, accZ_corr] = correctAccelerometer(imuTable, acc_offset)
    accX_corr = (imuTable.AccelX * 9.81 / 1000) - acc_offset(1);
    accY_corr = (imuTable.AccelY * 9.81 / 1000) - acc_offset(2);
    accZ_corr = (imuTable.AccelZ * 9.81 / 1000) - acc_offset(3);
    

end

function [GyroX_corr, GyroY_corr, GyroZ_corr] = correctGyroscope(imuTable, gyro_offset)
    GyroX_corr = deg2rad(imuTable.GyroX / 1000) - gyro_offset(1);
    GyroY_corr = deg2rad(imuTable.GyroY / 1000) - gyro_offset(2);
    GyroZ_corr = deg2rad(imuTable.GyroZ / 1000) - gyro_offset(3);

end