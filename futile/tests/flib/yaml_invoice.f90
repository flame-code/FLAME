!> @file
!! Test yaml output module
!! @example yaml_invoice.f90
!! Extended tests about yaml output generation
!! @author
!!    Copyright (C) 2013-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Test yaml part of flib
subroutine yaml_invoice_example()
  use yaml_output
  use yaml_strings, only: yaml_date_toa
  implicit none

  !call yaml_set_stream(tabbing=0)
  call yaml_comment('Yaml Invoice Example',hfill='-')
  call yaml_map('invoice',34843)
  call yaml_map('date',trim(yaml_date_toa()))
  call yaml_mapping_open('bill-to',label='id001')
   call yaml_map('given','Chris')
   call yaml_mapping_open('address')
      call yaml_mapping_open('lines')
      call yaml_scalar('458 Walkman Dr.')
      call yaml_scalar('Suite #292')
      call yaml_mapping_close()
   call yaml_mapping_close()
  call yaml_mapping_close()
  call yaml_map('ship_to','*id001')
  
  !next step: sequence elements
  call yaml_sequence_open('product')
  !call yaml_sequence_open()
    call yaml_sequence(advance='no')
!    call yaml_mapping_open()
      call yaml_map('sku','BL394D')
      call yaml_map('quantity',4)
      call yaml_map('description','Basketball')
      call yaml_map('price',450.,fmt='(f6.2)')
!    call yaml_mapping_close()
    !call yaml_newline() !new line in a flow 
     call yaml_sequence(advance='no')
!     call yaml_mapping_open()
     call yaml_map('sku','BL4438H')
     call yaml_map('quantity',1)
     call yaml_map('description','Super Hoop')
     call yaml_map('price',2392.,fmt='(f8.2)')
!     call yaml_mapping_close()
    call yaml_sequence_close()
    !final part
    call yaml_map('tax',251.42,fmt='(f6.2)')
    call yaml_map('total',4443.52d0,fmt='(f6.2)') !wrong format on purpose
    call yaml_map('comments','Late afternoon is best. Backup contact is Nancy Billsmer @ 338-4338.')

      !call yaml_mapping_close()

end subroutine yaml_invoice_example


!> Test yaml and dictionairies
subroutine yaml_invoice_example_with_dictionaries()
  use yaml_output
  use dictionaries
  use yaml_strings, only: yaml_date_toa
  implicit none
  real(kind=8) :: price
  type(dictionary), pointer :: dict,dict_tmp

  call yaml_comment('Yaml Invoice Example, using dictionaries',hfill='-')
  !setting the data in the fortran dictionary
  call dict_init(dict)
  call set(dict//'invoice',34843)
  call set(dict//'date',trim(yaml_date_toa()))

  call set(dict//'bill-to'//'given','Chris')
  call set(dict//'bill-to'//'family','Dumars')
  
  call dict_init(dict_tmp)
    call set(dict_tmp//'lines','458 Walkman Dr. Suite #292')
    call set(dict_tmp//'city','Royal Oak')
    call set(dict_tmp//'state','MI')
    call set(dict_tmp//'postal',48046)

  call set(dict//'bill-to'//'address',dict_tmp)
  !no need to free the dictionary after association
  
  !tagging of dictionary not yet implemented

    !products
  call dict_init(dict_tmp)
  call set(dict_tmp//'sku','BL34D')
  call set(dict_tmp//'quantity',4)
  call set(dict_tmp//'description','Basketball')
  call set(dict_tmp//'price',450.00)
  !adding to the item
  call add(dict//'Product',dict_tmp)

  call dict_init(dict_tmp)
  call set(dict_tmp//'sku','BL4438')
  call set(dict_tmp//'quantity',1)
  call set(dict_tmp//'description','Super Hoop')
  call set(dict_tmp//'price',2392.00)
  price=dict_tmp//'price'
  call yaml_map('Retrieve the price value',price)
  call add(dict//'Product',dict_tmp)

  call set(dict//'Tax',251.42)
  call set(dict//'Total',4443.52)
  call set(dict//'Comments','Late afternoon is best. Backup contact is Nancy Billsmer @ 338-4338')

  !print invoice
  call yaml_dict_dump(dict)

  call dict_free(dict)

end subroutine yaml_invoice_example_with_dictionaries
