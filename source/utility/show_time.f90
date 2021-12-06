!show date hour minute second
subroutine show_time()

integer, dimension(8)::time

call date_and_time(values=time)
write(*,'(I4, 1x, A4, 1x, I2, 1x, A5, 1x, I2, 1x, A3, 1x, I2, A1, I2, A1, I2)') &
time(1), 'year', time(2), 'month', time(3), 'day', time(5), ':', time(6), ':', time(7)

end subroutine show_time