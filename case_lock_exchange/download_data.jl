using HTTP
# download the locxx_map.nc data if needed
if ~isfile("locxx_map.nc")
    HTTP.download("https://nx7384.your-storageshare.de/s/fe2KHoLBGQx3S46/download", "locxx_map.nc")
end